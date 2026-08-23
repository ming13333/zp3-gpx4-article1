#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
#37 改进②：GPX4–ZP3 表达关联（回扣 Cell 2026 GPX4-ZP3 免疫抑制轴迁移逻辑）
=====================================================================
背景：三篇文章均声称将 Cell 2026「胞外 GPX4–ZP3 免疫抑制 DAMP 轴」从
生殖生物学迁移至神经肿瘤领域，但此前从未验证过 GPX4 与 ZP3 的表达关联。
本脚本用真实数据补齐这一回扣证据。

数据：本地真实 TPM (TcgaTargetGtex_rsem_gene_tpm.gz, 1.3GB, 已下载)
目标基因：
  - ZP3 (ENSG00000188372)
  - GPX4 (glutathione peroxidase 4, ENSG00000112715)  ← 轴核心
  - 铁死亡评分基因：GPX4 / ACSL4 / SLC7A11 / TFRC / FTL / FTH1 / NFE2L2 (NRF2)
    从 anti-ferroptosis(防御) 与 pro-ferroptosis(敏感) 两方面构成评分

分析内容：
  1. GPX4–ZP3 表达 Spearman 关联（全部 TCGA 肿瘤样本 + 胶质瘤 GBM/LGG）
  2. 铁死亡防御评分 (GPX4/SLC7A11/NFE2L2) 与 ZP3 关联
  3. 铁死亡敏感评分 (ACSL4/TFRC) 与 ZP3 关联
  4. 泛癌各癌种 GPX4–ZP3 关联分层，验证组织特异性
  5. 产出关联图 + CSV，回扣迁移逻辑

产物：gpx4_zp3_*csv + fig_gpx4_zp3_association.png
"""
import os, json, gzip, time
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(BASE)))
DATA = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling",
                    "TcgaTargetGtex_rsem_gene_tpm.gz")
DISEASE_MAP = os.path.join(ROOT, "output", "tcga_pancan", "tcga_disease_map.json")
ENSG_MAP = os.path.join(ROOT, "output", "tcga_pancan", "ensg_map.json")

# ---------------------------------------------------------------------------
# 目标基因（symbol -> ensg）
# ---------------------------------------------------------------------------
GENES = {
    "ZP3": "ENSG00000188372",
    # 铁死亡轴
    "GPX4": "ENSG00000112715",
    "ACSL4": "ENSG00000068366",
    "SLC7A11": "ENSG00000151012",
    "TFRC": "ENSG00000072274",      # 铁摄取促铁死亡
    "FTL": "ENSG00000087086",       # 铁蛋白轻链(储铁, 防铁死亡)
    "FTH1": "ENSG00000167996",      # 铁蛋白重链
    "NFE2L2": "ENSG00000116044",    # NRF2, 铁死亡防御转录因子
    # 下游脂质过氧化
    "ALOX15": "ENSG00000161905",
}

# 铁死亡评分分解
ANTI_FERROPTOSIS = ["GPX4", "SLC7A11", "NFE2L2", "FTL"]   # 防御
PRO_FERROPTOSIS = ["ACSL4", "TFRC"]                       # 敏感/促死


def ensure_ensg(dict_path, extra):
    """把缺失的基因映射合并写入缓存。"""
    with open(dict_path) as f:
        m = json.load(f)
    added = set()
    for sym, ensg in extra.items():
        if sym not in m:
            m[sym] = ensg
            added.add(sym)
    if added:
        with open(dict_path, "w") as f:
            json.dump(m, f, indent=2)
        print(f"  ensg_map.json 已补写 {sorted(added)}")
    return m


def read_target_genes(path, want_ensg, chunk=2000):
    """流式读取 gz TPM，仅保留目标基因行。返回 DataFrame: index=ensg, cols=sample。"""
    stripped = {e.split(".")[0]: e for e in want_ensg}
    rows = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        n = 0
        while True:
            lines = f.readlines(chunk)
            if not lines:
                break
            for ln in lines:
                parts = ln.rstrip("\n").split("\t")
                base = parts[0].split(".")[0]
                if base in stripped:
                    rows[stripped[base]] = [float(x) for x in parts[1:]]
            n += len(lines)
    return pd.DataFrame.from_dict(rows, orient="index", columns=samples)


def spearman_p(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 15:
        return np.nan, np.nan, int(m.sum())
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), int(m.sum())


def main():
    print("=== #37 GPX4–ZP3 表达关联分析（真实 TPM 数据）===\n")
    if not os.path.exists(DATA):
        print(f"!! 缺 TPM: {DATA}"); return

    # 1. 补写 ensg_map
    print("1. 确保基因映射（补 GPX4/ACSL4 等铁死亡基因）...")
    ensg_map = ensure_ensg(ENSG_MAP, GENES)
    want_ensg = list(GENES.values())
    inv = {v: k for k, v in GENES.items()}

    # 2. 读取目标基因
    print("2. 流式读取真实 TPM（目标基因）...")
    t0 = time.time()
    mat = read_target_genes(DATA, want_ensg)
    print(f"   矩阵 {mat.shape[0]} 基因 × {mat.shape[1]} 样本，耗时 {time.time()-t0:.1f}s")
    missing = [sym for sym, e in GENES.items() if e not in mat.index]
    if missing:
        print(f"   !! 缺失基因行: {missing}")
    mat.index = [inv.get(e, e) for e in mat.index]

    # 3. TCGA 肿瘤样本
    samples = list(mat.columns)
    tcga = [s for s in samples if s.startswith("TCGA-") and s.split("-")[3].startswith("01")]
    print(f"   TCGA 肿瘤样本: {len(tcga)}")
    mat_t = mat[tcga]

    # 4. 癌种映射（mat_t 是 基因行×样本列，转置后加 cancer 列）
    with open(DISEASE_MAP) as f:
        disease = json.load(f)
    def cancer_of(s):
        return disease.get("-".join(s.split("-")[:3]), "UNKNOWN")
    mat_tt = mat_t.T.copy()                # 样本行 × 基因列
    mat_tt["cancer"] = [cancer_of(s) for s in mat_tt.index]
    mat_tt = mat_tt[mat_tt["cancer"] != "UNKNOWN"]

    zp3 = mat_tt["ZP3"]
    gpx4 = mat_tt["GPX4"]

    # ---- 5. 全局 + 胶质瘤 GPX4–ZP3 ----
    print("\n3. GPX4–ZP3 关联：")
    rec_global = []
    for label, sel in [
        ("ALL_TCGA", np.ones(len(mat_tt), dtype=bool)),
        ("GBM", mat_tt["cancer"].values == "GBM"),
        ("LGG", mat_tt["cancer"].values == "LGG"),
    ]:
        rho, p, n = spearman_p(zp3.values[sel], gpx4.values[sel])
        rec_global.append({"Cohort": label, "Rho": rho, "P": p, "N": n})
        print(f"   {label:10s}: GPX4–ZP3 ρ={rho:+.3f}, p={p:.3g}, n={n}")

    # ---- 6. 铁死亡评分与 ZP3 ----
    def score_zs(genes, df):
        """z-score 共识：每基因（列）跨样本标准化后，对每个样本取基因均值。"""
        sub = df[genes].astype(float).T        # 基因行 × 样本列
        valid = sub.std(axis=1) > 0            # 过滤 std=0 基因
        if not valid.any():
            return pd.Series(np.nan, index=df.index)
        z = ((sub.loc[valid] - sub.loc[valid].mean(axis=1).values[:, None])
             / sub.loc[valid].std(axis=1).values[:, None])
        return z.mean(axis=0)                  # 每样本一个评分

    rec_fx = []
    for label, gs in [("AntiFerroptosis", ANTI_FERROPTOSIS),
                      ("ProFerroptosis", PRO_FERROPTOSIS)]:
        available = [g for g in gs if g in mat_tt.columns]
        if not available:
            continue
        s = score_zs(available, mat_tt)
        rho, p, n = spearman_p(zp3.values, s.values)
        rec_fx.append({"Score": label, "Genes": "+".join(available),
                       "Rho": rho, "P": p, "N": n})
        print(f"   {label:16s}: ZP3–{label} ρ={rho:+.3f}, p={p:.3g}, n={n}")

    # ---- 7. 泛癌各癌种 GPX4–ZP3 分层 ----
    rec_pan = []
    for c, grp in mat_tt.groupby("cancer"):
        if len(grp) < 30:
            continue
        rho, p, n = spearman_p(grp["ZP3"].values, grp["GPX4"].values)
        rec_pan.append({"Cancer": c, "Rho": rho, "P": p, "N": n})
    pan = pd.DataFrame(rec_pan).sort_values("Rho", ascending=False)
    print("\n4. 泛癌 GPX4–ZP3 关联（Top 10）:")
    for _, r in pan.head(10).iterrows():
        star = "*" if r["P"] < 0.05 else " "
        print(f"   {r['Cancer']:6s}  ρ={r['Rho']:+.3f}  p={r['P']:.3g}  n={r['N']} {star}")

    # ---- 8. 保存 CSV ----
    pd.DataFrame(rec_global).to_csv(os.path.join(BASE, "gpx4_zp3_global.csv"), index=False)
    pd.DataFrame(rec_fx).to_csv(os.path.join(BASE, "gpx4_zp3_ferroptosis_score.csv"), index=False)
    pan.to_csv(os.path.join(BASE, "gpx4_zp3_pancancer.csv"), index=False)
    mat_tt[["ZP3", "GPX4", "cancer"]].to_csv(
        os.path.join(BASE, "gpx4_zp3_expr_matrix.csv"))

    # ---- 9. 图 ----
    print("\n5. 生成关联图...")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))

    # a) 全 TCGA 散点
    ax = axes[0]
    ax.scatter(gpx4.values, zp3.values, s=8, alpha=0.4, c="#378ADD")
    rho_all, p_all, _ = spearman_p(zp3.values, gpx4.values)
    ax.set_xlabel("GPX4 log2(TPM)")
    ax.set_ylabel("ZP3 log2(TPM)")
    ax.set_title(f"All TCGA (n={len(mat_tt)})\nρ={rho_all:+.3f}, p={p_all:.1e}", fontsize=11)

    # b) GBM/LGG 散点
    ax = axes[1]
    for c, col in [("GBM", "#A32D2D"), ("LGG", "#1D9E75")]:
        g = mat_tt[mat_tt["cancer"] == c]
        ax.scatter(g["GPX4"].values, g["ZP3"].values, s=12, alpha=0.6,
                   label=c, c=col)
    ax.set_xlabel("GPX4 log2(TPM)"); ax.set_ylabel("ZP3 log2(TPM)")
    ax.set_title("Glioma: GPX4–ZP3", fontsize=11)
    ax.legend()

    # c) 泛癌条形
    ax = axes[2]
    top = pan.head(12)
    colors = ["#A32D2D" if (r["Rho"] > 0 and r["P"] < 0.05) else "#888780"
              for _, r in top.iterrows()]
    ax.barh(range(len(top)), top["Rho"], color=colors)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["Cancer"])
    ax.axvline(0, color="gray", lw=0.8)
    ax.set_xlabel("GPX4–ZP3 Spearman ρ")
    ax.set_title("GPX4–ZP3 by cancer type", fontsize=11)

    plt.tight_layout()
    fig.savefig(os.path.join(BASE, "fig_gpx4_zp3_association.png"),
                dpi=300, bbox_inches="tight")
    print("   已保存 fig_gpx4_zp3_association.png")

    # ---- 10. 回扣迁移逻辑结论 ----
    print("\n=== 结论（回扣 Cell 2026 GPX4-ZP3 轴迁移逻辑）===")
    gbm_r = next(r for r in rec_global if r["Cohort"] == "GBM")
    lgg_r = next(r for r in rec_global if r["Cohort"] == "LGG")
    all_r = next(r for r in rec_global if r["Cohort"] == "ALL_TCGA")
    print(f"  · 全局 TCGA: GPX4–ZP3 ρ={all_r['Rho']:+.3f} (p={all_r['P']:.2g})")
    print(f"  · GBM: ρ={gbm_r['Rho']:+.3f} (p={gbm_r['P']:.2g})")
    print(f"  · LGG: ρ={lgg_r['Rho']:+.3f} (p={lgg_r['P']:.2g})")
    for r in rec_fx:
        print(f"  · ZP3–{r['Score']}: ρ={r['Rho']:+.3f} (p={r['P']:.2g})")
    print("\n回扣结论：GPX4 作为铁死亡轴核心，若与 ZP3 显著正相关，"
          "则支持『铁死亡(ferroptosis) 释放 GPX4→分泌型 GPX4-ZP3 复合物→"
          "髓系免疫抑制』这一 Cell 2026 轴在神经肿瘤中的延伸。"
          "若 GBM/LGG 与泛癌方向一致，则为迁移逻辑提供表达层面证据。")
    print("\n=== #37 完成 ===")


if __name__ == "__main__":
    main()
