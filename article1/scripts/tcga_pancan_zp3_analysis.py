#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TCGA 泛癌分析 ZP3 与免疫特征的组织特异性（真实数据版）
=================================================================
数据：UCSC Xena「TCGA TARGET GTEx」RSEM gene TPM
      (TcgaTargetGtex_rsem_gene_tpm.gz，行=Ensembl 基因 ID，列=样本)
      —— 该文件由 segmented_resume_download 真实下载，非模拟数据。

流程：
  1. 解析 ZP3 与各免疫标志基因 symbol -> Ensembl ID（Ensembl REST 获取并缓存到
     ensg_map.json，离线可复现；失败则用内置回退表）。
  2. 流式读取 gz 文件，仅保留目标基因行（避免把 1.3GB 全量读入内存）。
  3. 筛选 TCGA 肿瘤样本（sample id 形如 TCGA-XXXX-####-01），按癌种分组。
  4. 对每个癌种、每个免疫基因集，计算 ZP3 表达与「基因集均值评分」的
     Spearman 秩相关（表达量为 log2(TPM) 右侧偏态，用 Spearman 而非 Pearson）。
  5. 生成泛癌热图 + 各癌种关联强度柱状图，并保存明细/汇总 CSV。

注：表达值为 log2(TPM)；"-9.9658" 为 Xena 对 0 的下溢底值，保留为低表达秩
    （秩相关对单调变换稳健，不参与绝对量解释）。
"""
import os, sys, json, gzip, time
from urllib.parse import urlencode
import urllib.request
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(BASE)))  # 项目根
DATA = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling",
                    "TcgaTargetGtex_rsem_gene_tpm.gz")
CACHE = os.path.join(BASE, "ensg_map.json")

# ---------------------------------------------------------------------------
# 1. 免疫基因集（symbol 定义，与旧版一致）
# ---------------------------------------------------------------------------
IMMUNE_GENE_SETS = {
    'M2_Macrophage': ['CD163', 'MSR1', 'MRC1', 'VSIG4', 'CD200R1', 'TGFB1', 'IL10',
                      'ARG1', 'MERTK', 'CLEC7A'],
    'T_cell_exhaustion': ['LAG3', 'TIGIT', 'HAVCR2', 'PDCD1', 'CTLA4', 'CD274',
                          'PDCD1LG2', 'BTLA', 'VSIR', 'IDO1', 'IDO2'],
    'Cytolytic_activity': ['GZMA', 'GZMB', 'PRF1', 'IFNG'],
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF18', 'ICOS', 'CD40LG'],
    'IFN_gamma': ['IFNG', 'STAT1', 'IRF1', 'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'CD274'],
    'Checkpoint': ['CD274', 'PDCD1', 'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2', 'BTLA', 'VSIR'],
    'Myeloid': ['CD68', 'CD163', 'CSF1R', 'ITGAM', 'CD14', 'LYZ', 'S100A8', 'S100A9'],
}
ZP3_SYMBOL = "ZP3"

# 内置回退表（仅当 Ensembl REST 不可用时使用，覆盖最关键基因）
FALLBACK_ENSG = {
    "ZP3": "ENSG00000188372", "CD163": "ENSG00000177697", "CD68": "ENSG00000097273",
    "CD274": "ENSG00000120217", "PDCD1": "ENSG00000188389", "CTLA4": "ENSG00000163558",
    "FOXP3": "ENSG00000049768", "IFNG": "ENSG00000111537", "TGFB1": "ENSG00000105329",
    "IL10": "ENSG00000136634", "TIGIT": "ENSG00000180554", "LAG3": "ENSG00000179868",
    "VSIR": "ENSG00000116286", "BTLA": "ENSG00000120735", "IDO1": "ENSG00000131203",
    "TREM2": "ENSG00000183463",
}


def get_ensg_map(symbols):
    """symbol -> Ensembl gene id（去版本号）。优先读缓存；否则查 Ensembl REST 并缓存。"""
    if os.path.exists(CACHE):
        with open(CACHE) as f:
            cache = json.load(f)
    else:
        cache = {}
    needed = [s for s in symbols if s not in cache]
    if needed:
        import urllib.request
        print(f"  从 Ensembl REST 解析 {len(needed)} 个 symbol -> ENSG ...")
        for s in needed:
            try:
                url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{s}?content-type=application/json"
                req = urllib.request.Request(url, headers={"User-Agent": "wb-research/1.0"})
                with urllib.request.urlopen(req, timeout=20) as r:
                    j = json.loads(r.read().decode())
                if "id" in j:
                    cache[s] = j["id"].split(".")[0]
            except Exception as e:
                fb = FALLBACK_ENSG.get(s)
                if fb:
                    cache[s] = fb
                else:
                    print(f"    !! {s} 无法解析: {e}")
            time.sleep(0.05)
        with open(CACHE, "w") as f:
            json.dump(cache, f, indent=2)
    return {s: cache.get(s) for s in symbols}


# ---------------------------------------------------------------------------
# 2. 流式读取目标基因
# ---------------------------------------------------------------------------
def get_tcga_disease_map():
    """TCGA 参与者 barcode -> 癌种（如 BRCA/LUAD/LGG/GBM）。
    样本 ID 形如 TCGA-<TSS>-<参与者>-<样本>，其中第 2 段是组织来源地(TSS)，
    并非癌种；癌种需从 GDC 的 case->project 映射获得。结果缓存到
    tcga_disease_map.json，离线可复现。"""
    cache = os.path.join(BASE, "tcga_disease_map.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return json.load(f)
    cases = {}
    base = "https://api.gdc.cancer.gov/cases"
    size, frm = 1000, 0
    while True:
        q = urlencode({
            "fields": "submitter_id,project.project_id",
            "size": size, "from": frm})
        try:
            req = urllib.request.Request(base + "?" + q,
                                         headers={"User-Agent": "wb-research/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode())
        except Exception as e:
            print(f"  !! GDC 映射获取失败: {e}")
            break
        hits = data.get("data", {}).get("hits", [])
        if not hits:
            break
        for h in hits:
            sid = h.get("submitter_id")
            pid = (h.get("project") or {}).get("project_id")
            if sid and pid and pid.startswith("TCGA-"):
                cases[sid] = pid.split("-")[-1]  # TCGA-BRCA -> BRCA
        frm += size
        if len(hits) < size:
            break
        time.sleep(0.05)  # 避免触发 GDC 限流
    with open(cache, "w") as f:
        json.dump(cases, f, indent=2)
    print(f"  GDC 映射缓存 {len(cases)} 个参与者->癌种")
    return cases


def read_target_genes(path, target_ids, chunk=2000):
    """流式读取 gz TPM 矩阵，仅保留 target_ids（Ensembl，含或不含版本）。
    返回 DataFrame：index=Ensembl id，columns=sample id。"""
    target_strip = {t.split(".")[0]: t for t in target_ids}
    rows = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        kept_samples = samples  # 先全保留，后面再筛 TCGA
        n = 0
        while True:
            lines = f.readlines(chunk)
            if not lines:
                break
            for ln in lines:
                parts = ln.rstrip("\n").split("\t")
                gid = parts[0]
                base = gid.split(".")[0]
                if base in target_strip:
                    rows[target_strip[base]] = [float(x) for x in parts[1:]]
            n += len(lines)
            if n % 10000 == 0:
                print(f"    已扫描 {n} 个基因行...")
    df = pd.DataFrame.from_dict(rows, orient="index", columns=kept_samples)
    return df


# ---------------------------------------------------------------------------
# 3. 主流程
# ---------------------------------------------------------------------------
def main():
    print("=== TCGA 泛癌分析 ZP3 与免疫特征的组织特异性（真实数据）===\n")
    if not os.path.exists(DATA):
        print(f"!! 找不到真实数据文件: {DATA}\n   请先运行 segmented_resume_download 下载 "
              f"TcgaTargetGtex_rsem_gene_tpm.gz")
        sys.exit(1)

    all_symbols = [ZP3_SYMBOL] + sorted({g for s in IMMUNE_GENE_SETS.values() for g in s})
    sym2ensg = get_ensg_map(all_symbols)
    unresolved = [s for s, e in sym2ensg.items() if not e]
    if unresolved:
        print(f"  !! 未解析的基因（将跳过）: {unresolved}")
    # 目标 Ensembl id 列表
    target_ids = [e for e in sym2ensg.values() if e]
    print(f"  解析到 {len(target_ids)} 个 Ensembl 基因（含 ZP3）\n")

    print("2. 流式读取真实 TPM 矩阵（仅目标基因）...")
    t0 = time.time()
    mat = read_target_genes(DATA, target_ids)
    print(f"  读取完成，矩阵 {mat.shape[0]} 基因 × {mat.shape[1]} 样本，"
          f"耗时 {time.time()-t0:.1f}s")

    # ZP3 行
    zp3_ensg = sym2ensg.get(ZP3_SYMBOL)
    if zp3_ensg not in mat.index:
        print(f"!! 矩阵中无 ZP3 ({zp3_ensg})，退出"); sys.exit(1)
    zp3_vec_all = mat.loc[zp3_ensg]

    # 筛选 TCGA 肿瘤样本
    samples = list(mat.columns)
    tcga_mask = [s.startswith("TCGA-") and s.split("-")[3].startswith("01") for s in samples]
    tcga_samples = [s for s, m in zip(samples, tcga_mask) if m]
    print(f"  TCGA 肿瘤样本数: {len(tcga_samples)}")
    mat_t = mat[tcga_samples]

    # 按癌种分组：样本 barcode = TCGA-<TSS>-<参与者>-<样本>，
    # 第2段是组织来源地(TSS)，非癌种；用 GDC 参与者->癌种映射。
    participant_of = {s: "-".join(s.split("-")[:3]) for s in tcga_samples}
    p2cancer = get_tcga_disease_map()  # 参与者 barcode -> 癌种缩写
    cancer_of = {}
    for s in tcga_samples:
        cancer_of[s] = p2cancer.get(participant_of[s], "UNKNOWN")
    cancers = {}
    for s in tcga_samples:
        cancers.setdefault(cancer_of[s], []).append(s)
    cancers = {c: v for c, v in cancers.items() if c != "UNKNOWN"}
    n_unknown = sum(1 for c in cancer_of.values() if c == "UNKNOWN")
    if n_unknown:
        print(f"  (警告：{n_unknown} 个样本未能映射癌种，已排除)")

    print(f"\n3. 计算各癌种 ZP3-免疫关联（Spearman rho）...")
    print("   评分方法：z-score 共识法（每基因跨样本标准化后取均值）")
    records = []
    for cancer, sams in cancers.items():
        if len(sams) < 30:
            continue
        zp3_v = mat_t.loc[zp3_ensg, sams].values.astype(float)
        for set_name, syms in IMMUNE_GENE_SETS.items():
            ensgs = [sym2ensg[s] for s in syms if sym2ensg.get(s) in mat_t.index]
            if not ensgs:
                continue
            # z-score 共识法：每基因先跨样本标准化，再取均值
            sub = mat_t.loc[ensgs, sams]           # 基因×样本
            gene_mean = sub.mean(axis=1)           # 每基因均值
            gene_std = sub.std(axis=1)             # 每基因标准差
            # 过滤 std=0 的基因（如全为底值）
            valid = gene_std > 0
            if not valid.any():
                continue
            z = (sub.loc[valid] - gene_mean[valid].values[:, None]) / gene_std[valid].values[:, None]
            score = z.mean(axis=0).values.astype(float)
            mask = np.isfinite(zp3_v) & np.isfinite(score)
            if mask.sum() < 20:
                continue
            rho, p = stats.spearmanr(zp3_v[mask], score[mask])
            records.append({
                "Cancer_Code": cancer, "Feature": set_name,
                "Rho": round(float(rho), 4), "P_value": float(p),
                "N": int(mask.sum()),
                "Significant": bool(p < 0.05),
            })
    res = pd.DataFrame(records)
    if res.empty:
        print("!! 无有效关联结果（检查基因解析/样本量）"); sys.exit(1)
    print(f"  共 {len(res)} 条 (癌种×基因集) 关联结果，覆盖 {res['Cancer_Code'].nunique()} 种癌")

    # 各癌种平均关联
    cancer_summary = res.groupby("Cancer_Code").agg(
        Avg_Rho=("Rho", "mean"),
        Sig_Count=("Significant", "sum"),
        N_sets=("Feature", "count"),
    ).reset_index().sort_values("Avg_Rho", ascending=False)
    cancer_summary["Sig_Count"] = cancer_summary["Sig_Count"].astype(int)

    print("\n  ZP3-免疫关联强度排名 (Top 8):")
    print("  " + "-" * 58)
    for _, row in cancer_summary.head(8).iterrows():
        print(f"  {row['Cancer_Code']:6s} | ρ = {row['Avg_Rho']:+.3f} | "
              f"Sig sets: {row['Sig_Count']}/{int(row['N_sets'])} | n={row['N_sets']*0}")

    # 热图
    print("\n4. 生成泛癌关联热图与柱状图...")
    heat = res.pivot_table(index="Cancer_Code", columns="Feature", values="Rho", aggfunc="mean")
    heat = heat.loc[cancer_summary["Cancer_Code"]]
    fig, axes = plt.subplots(1, 2, figsize=(18, 12))
    sns.heatmap(heat, cmap="RdBu_r", center=0, vmin=-0.3, vmax=0.3,
                ax=axes[0], cbar_kws={"label": "Spearman ρ"})
    axes[0].set_title("ZP3–Immune Feature Correlations\nAcross TCGA Cancer Types (real data)")
    axes[0].set_xlabel("Immune Feature (gene-set score)")
    axes[0].set_ylabel("Cancer Type")

    top = cancer_summary.head(15)
    colors = ["#e74c3c" if r > 0.15 else "#f39c12" if r > 0.05 else "#3498db"
              for r in top["Avg_Rho"]]
    axes[1].barh(range(len(top)), top["Avg_Rho"], color=colors)
    axes[1].set_yticks(range(len(top)))
    axes[1].set_yticklabels(top["Cancer_Code"])
    axes[1].axvline(0.05, color="gray", ls="--", alpha=0.5, label="weak threshold")
    axes[1].axvline(0.15, color="gray", ls="-", alpha=0.5, label="moderate threshold")
    axes[1].set_xlabel("Mean Spearman ρ")
    axes[1].set_title("ZP3–Immune Association Strength\nby Cancer Type (Top 15, real data)")
    axes[1].legend()
    plt.tight_layout()
    fig.savefig(os.path.join(BASE, "fig_tcga_pancan_zp3_heatmap.png"), dpi=300, bbox_inches="tight")
    print("  已保存 fig_tcga_pancan_zp3_heatmap.png")

    res.to_csv(os.path.join(BASE, "tcga_pancan_zp3_correlations.csv"), index=False)
    cancer_summary.to_csv(os.path.join(BASE, "tcga_pancan_cancer_summary.csv"), index=False)
    print("  已保存 tcga_pancan_zp3_correlations.csv / tcga_pancan_cancer_summary.csv")

    # 结论
    print("\n5. 结论（真实数据）:")
    strong = cancer_summary[cancer_summary["Avg_Rho"] > 0.15]["Cancer_Code"].tolist()
    moderate = cancer_summary[(cancer_summary["Avg_Rho"] > 0.05) &
                              (cancer_summary["Avg_Rho"] <= 0.15)]["Cancer_Code"].tolist()
    weak = cancer_summary[cancer_summary["Avg_Rho"] <= 0.05]["Cancer_Code"].tolist()
    print(f"  强关联 (ρ>0.15): {len(strong)} 种 — {', '.join(strong[:8])}")
    print(f"  中等关联 (0.05<ρ≤0.15): {len(moderate)} 种 — {', '.join(moderate[:8])}")
    print(f"  弱/无关联 (ρ≤0.05): {len(weak)} 种 — {', '.join(weak[:8])}")
    print("\n=== 分析完成（真实 TCGA TARGET GTEx 数据，已替代原模拟数据）===")


if __name__ == "__main__":
    main()
