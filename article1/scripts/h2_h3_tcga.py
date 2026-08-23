# -*- coding: utf-8 -*-
"""
H2 (预后关联) + H3 (免疫抑制 TME 关联) — TCGA GBM/LGG bulk RNA-seq
数据：UCSC Xena 下载的 HiSeq_TCGA_gene.xena.gz (表达, 行=基因, 列=样本)
      GBM/LGG_clinicalMatrix.gz (临床, 含 OS_MONTHS / OS_STATUS)
描述级分析，诚实标注。
"""
import os, gzip, sys, numpy as np, pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))

# 统一使用经 lifelines/标准实现验证的 log-rank（消除旧版两套不一致手写实现）
sys.path.insert(0, os.path.join(BASE, "..", "common"))
from stats_utils import logrank  # noqa: E402

# ---- 免疫抑制相关标志基因（文献常见，作描述性关联）----
IMMUNOSUPP_GENES = ["TGFB1", "IL10", "FOXP3", "CD274", "PDCD1", "CTLA4",
                    "MRC1", "CD163", "VSIG4", "ARG1", "IDO1", "VEGFA",
                    "CCL2", "CXCL12", "MSR1", "TREM2"]
# M2 / TAM 偏向
M2_GENES = ["MRC1", "CD163", "MSR1", "ARG1", "TGFB1", "IL10", "VSIG4"]
# Treg
TREG_GENES = ["FOXP3", "IL2RA", "CTLA4", "TIGIT"]
# 检查点
CHECKPT_GENES = ["CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3"]

def load_expr(path):
    """Xena 表达矩阵：行=基因符号, 列=样本。返回 DataFrame。"""
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    return df

def load_clin(path):
    """Xena 临床矩阵：尝试解析为 属性×样本 或 样本×属性。"""
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    # 若 OS_MONTHS 在 index -> 属性×样本
    if "OS_MONTHS" in df.index or "OS_MONTHS" in df.columns:
        if "OS_MONTHS" in df.index:
            return df.T  # -> 样本×属性
        return df
    return df

def extract_os(clin, sample_ids):
    """从临床表提取 OS_MONTHS / OS_STATUS，匹配表达样本。"""
    # clin: 样本×属性
    common = [s for s in sample_ids if s in clin.index]
    sub = clin.loc[common]
    # 找 OS 列（灵活匹配）
    os_time_col = next((c for c in sub.columns if c.upper() in ("OS_MONTHS", "_OS_MONTHS", "OS_MONTHS")), None)
    os_stat_col = next((c for c in sub.columns if "OS_STATUS" in c.upper() or c.upper() in ("_OS_STATUS",)), None)
    if os_time_col is None or os_stat_col is None:
        # 打印可用列帮助调试
        print("  临床可用列(前30):", list(sub.columns[:30]))
        return None
    t = pd.to_numeric(sub[os_time_col], errors="coerce")
    # OS_STATUS: 'DECEASED'/'1' = event=1; 'LIVING'/'0' = 0
    s = sub[os_stat_col].astype(str).str.upper()
    event = s.map(lambda x: 1 if ("DECEAS" in x or x.strip() == "1") else (0 if ("LIV" in x or x.strip() == "0") else np.nan))
    out = pd.DataFrame({"time": t, "event": event}, index=sub.index)
    out = out.dropna()
    return out

def analyze_cohort(name, expr_path, clin_path):
    print("\n" + "=" * 72)
    print(f"队列: {name}")
    expr = load_expr(expr_path)
    clin = load_clin(clin_path)
    print(f"  表达矩阵: {expr.shape[0]} 基因 × {expr.shape[1]} 样本")
    if "ZP3" not in expr.index:
        print("  !! 表达矩阵无 ZP3 行，跳过")
        return
    zp3 = expr.loc["ZP3"]
    # 只保留数值样本
    zp3 = pd.to_numeric(zp3, errors="coerce").dropna()
    osd = extract_os(clin, zp3.index.tolist())
    if osd is None or len(osd) < 20:
        print("  !! OS 数据不足，跳过 survival")
        return
    # 合并
    merged = pd.concat([zp3.rename("ZP3"), osd], axis=1).dropna()
    merged = merged[merged["time"] > 0]
    print(f"  survival 可用样本: {len(merged)}")
    # 中位二分
    med = merged["ZP3"].median()
    merged["group"] = (merged["ZP3"] > med).astype(int)
    hi = merged[merged.group == 1]; lo = merged[merged.group == 0]
    chi2, p = logrank(merged["time"].values, merged["event"].values, merged["group"].values)
    print(f"  ZP3 中位={med:.3f} | High n={len(hi)} vs Low n={len(lo)}")
    print(f"  H2 logrank: chi2={chi2:.3f}, p={p:.4g}")
    # 方向：High 组事件率
    rate_hi = hi["event"].mean(); rate_lo = lo["event"].mean()
    print(f"  事件率 High={rate_hi:.2f} Low={rate_lo:.2f} -> {'High ZP3 预后更差' if rate_hi>rate_lo else 'High ZP3 预后更好'}")

    # H3: ZP3 与免疫抑制标志关联
    # 注：表达量为右侧偏态、非正态，改用 Spearman 秩相关（原为 Pearson，已修正）
    print(f"  --- H3: ZP3 vs 免疫抑制标志 (Spearman rho) ---")
    rows = []
    h3_genes = list(dict.fromkeys(IMMUNOSUPP_GENES + M2_GENES + TREG_GENES + CHECKPT_GENES))  # 去重（同一基因在多集合）
    for gene in h3_genes:
        if gene in expr.index:
            g = pd.to_numeric(expr.loc[gene], errors="coerce")
            gg = pd.concat([merged["ZP3"], g.rename(gene)], axis=1).dropna()
            if len(gg) > 20:
                r, pp = stats.spearmanr(gg["ZP3"], gg[gene])
                rows.append((gene, round(r, 3), round(pp, 4), len(gg)))
    h3 = pd.DataFrame(rows, columns=["gene", "spearman_rho", "p", "n"]).sort_values("spearman_rho", ascending=False)
    print(h3.to_string(index=False))
    # 复合免疫抑制评分（z-mean of M2+TREG+CHECKPT available genes）
    sig = [g for g in M2_GENES + TREG_GENES + CHECKPT_GENES if g in expr.index]
    if sig:
        sub = expr.loc[sig, merged.index].apply(pd.to_numeric, errors="coerce")
        z = (sub - sub.mean()) / sub.std()
        immuno_score = z.mean(axis=0)
        gg = pd.concat([merged["ZP3"], immuno_score.rename("immuno_score")], axis=1).dropna()
        r, pp = stats.spearmanr(gg["ZP3"], gg["immuno_score"])
        print(f"  复合免疫抑制评分 vs ZP3: rho={r:.3f}, p={pp:.4g}, n={len(gg)}")
    return merged, h3

if __name__ == "__main__":
    results = {}
    results["GBM"] = analyze_cohort("TCGA GBM",
        os.path.join(BASE, "TCGA.GBM.sampleMap/HiSeq_TCGA_gene.xena.gz"),
        os.path.join(BASE, "TCGA.GBM.sampleMap/GBM_clinicalMatrix.gz"))
    results["LGG"] = analyze_cohort("TCGA LGG",
        os.path.join(BASE, "TCGA.LGG.sampleMap/HiSeq_TCGA_gene.xena.gz"),
        os.path.join(BASE, "TCGA.LGG.sampleMap/LGG_clinicalMatrix.gz"))
    print("\n=== H2/H3 分析完成（描述级，需外部验证）===")
