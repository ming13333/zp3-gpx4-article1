#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Article 1 补充图：ZP3 亚型方向链（跨队列一致性）
===================================================
核心亮点：#38 发现"亚型方向链"——ZP3 沿预后恶化方向递增：
    G-CIMP(预后最好/IDH-Mut) → Codel(1p19q共缺失) → Non-codel → IDH-WT(预后最差)

本图整合 3 个独立队列：
  · TCGA-GBM     Verhaak 亚型（G-CIMP 组 = IDH-Mut，预后最好）
  · TCGA-LGG     1p19q Codel/Non-codel + IDH 状态
  · CGGA-693     IDH_mutation_status + 1p19q_codeletion_status（独立队列）
                  以及 IDH 状态（Wildtype/Mutant）

为跨队列可比，Y 轴采用"相对本队列最低组的倍数变化"（FPKM/TPM 单位不同不可直接比绝对值），
并在次级标注各队列的绝对 ZP3 中位数。

产物：
  tcga_subtype_direction_chain.png
  subtype_direction_chain_values.csv
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
CGG = os.path.join(os.path.dirname(BASE), "cgga_validation")  # output/cgga_validation
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article1", "results", "tcga_subtype_results")
os.makedirs(OUT, exist_ok=True)

# ---------- 1. 载入 #38 已产出的 TCGA 亚型数据 ----------
TCGA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article1", "results", "tcga_subtype_zp3_results")
gbm = pd.read_csv(os.path.join(TCGA_DIR, "gbm_zp3_subtype_annotated.csv"))
lgg = pd.read_csv(os.path.join(TCGA_DIR, "lgg_zp3_subtype_annotated.csv"))

# ---------- 2. 载入 CGGA-693 ----------
cl = pd.read_csv(os.path.join(CGG, "CGGA.mRNAseq_693_clinical.20200506.txt"), sep="\t")
cl.columns = [c.split()[0] for c in cl.columns]
expr = pd.read_csv(os.path.join(CGG, "CGGA.mRNAseq_693.RSEM-genes.20200506.txt"),
                   sep="\t", index_col=0)
zp3_cgga = expr.loc["ZP3"]
cgga_df = cl.set_index("CGGA_ID")[["IDH_mutation_status", "1p19q_codeletion_status"]].copy()
cgga_df["ZP3"] = zp3_cgga
cgga_df = cgga_df[cgga_df["ZP3"].notna()]

# ---------- 3. 组装方向链（每级 = (标签, 队列, 该级全部 ZP3 值)） ----------
# TCGA-GBM 亚型（G-CIMP 用中位数，注意单位 log2 TPM）
gbm_gcimp = gbm[gbm["subtype"] == "G-CIMP"]["ZP3"]
gbm_non_gcimp = gbm[gbm["subtype"] != "G-CIMP"]["ZP3"]

# TCGA-LGG
lgg_codel = lgg[lgg["idh_1p19q"] == "Codel"]["ZP3"]
lgg_noncodel = lgg[lgg["idh_1p19q"] == "Non-codel"]["ZP3"]
lgg_wt = lgg[lgg["idh_status"] == "WT"]["ZP3"]
lgg_mut = lgg[lgg["idh_status"] == "Mutant"]["ZP3"]

# CGGA-693（FPKM 单位）
cg_idh_wt = cgga_df[cgga_df["IDH_mutation_status"] == "Wildtype"]["ZP3"]
cg_idh_mut = cgga_df[cgga_df["IDH_mutation_status"] == "Mutant"]["ZP3"]
cg_codel = cgga_df[cgga_df["1p19q_codeletion_status"] == "Codel"]["ZP3"]
cg_noncodel = cgga_df[cgga_df["1p19q_codeletion_status"] == "Non-codel"]["ZP3"]

def chain(label, cohort, series, color):
    return {"label": label, "cohort": cohort, "values": series.dropna().values,
            "color": color}

# 使用"分子病理方向链"（非转录组亚型，避免 GBM/LGG 转录组方向矛盾）：
# 排序关键：IDH-Mut（预后好）在左，IDH-WT（预后差）在右
rows = [
    # IDH-Mut 组（低 ZP3）
    chain("GBM\nG-CIMP",   "TCGA", gbm_gcimp,   "#1D9E75"),
    chain("LGG\nCodel",    "TCGA", lgg_codel,   "#2E8B57"),
    chain("CGGA\nCodel",   "CGGA", cg_codel,    "#7BC88A"),
    chain("LGG\nNon-codel","TCGA", lgg_noncodel,"#888780"),
    chain("CGGA\nNon-Codel","CGGA",cg_noncodel, "#B8B8B8"),
    chain("LGG\nIDH-Mut",  "TCGA", lgg_mut,     "#378ADD"),
    # IDH-WT 组（高 ZP3）
    chain("CGGA\nIDH-WT",  "CGGA", cg_idh_wt,   "#A32D2D"),
    chain("LGG\nIDH-WT",   "TCGA", lgg_wt,      "#D9534F"),
]

# ---------- 4. 相对本队列最低组倍数变化 ----------
# 每个队列内：选中一组作为基线（TCGA 用 G-CIMP；CGGA 用 Codel 或 IDH-Mut 中较低者）
def cohort_baseline(cohort, rows):
    """返回该队列所有组中 ZP3 中位数最低者。"""
    vals = [np.median(r["values"]) for r in rows if r["cohort"] == cohort and len(r["values"])]
    return min(vals) if vals else 1.0

baselines = {c: cohort_baseline(c, rows) for c in ("TCGA", "CGGA")}

rec = []
for r in rows:
    if len(r["values"]) < 5:
        continue
    med = float(np.median(r["values"]))
    fc = med / baselines[r["cohort"]] if baselines[r["cohort"]] > 0 else 1.0
    rec.append({
        "group": r["label"].replace("\n", " "), "cohort": r["cohort"],
        "median_ZP3": med, "base_FC": fc, "N": len(r["values"]), "color": r["color"]})
val_df = pd.DataFrame(rec)
val_df.to_csv(os.path.join(OUT, "subtype_direction_chain_values.csv"), index=False)

# ---------- 5. 绘图 ----------
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, 5.2),
                              gridspec_kw={"width_ratios": [3, 1.2]})

# 主面板：方向链（按原顺序）—— X 用序数位置，按 IDH-Mut→WT 排列
x = np.arange(len(rec))
ax.bar(x, val_df["base_FC"], color=val_df["color"], alpha=0.85,
       edgecolor="white", linewidth=0.8)
ax.axhline(1.0, color="#666", lw=0.8, ls="--")
ax.set_xticks(x)
ax.set_xticklabels(val_df["group"], fontsize=9)
ax.set_ylabel("Median ZP3  /  cohort-min median\n(fold-change, log-scaled)")
ax.set_yscale("log")
ax.set_title("ZP3 'molecular-subtype direction chain' across cohorts\n"
             "IDH-mutated / low-grade (better prognosis) → IDH-WT (worse prognosis)",
             fontsize=11)
# 标注绝对中位数
for xi, row in zip(x, val_df.itertuples()):
    unit = "log2TPM" if row.cohort == "TCGA" else "FPKM"
    ax.text(xi, row.base_FC * 1.06, f"{row.median_ZP3:.2f}\n({unit})",
            ha="center", va="bottom", fontsize=7.5, color="#222")

# 右侧：IDH-Mutant vs Wildtype 跨队列倍数对比
idh_pairs = [
    ("TCGA", lgg_mut, lgg_wt, "LGG-IDH"),
    ("CGGA", cg_idh_mut, cg_idh_wt, "CGGA-IDH"),
]
fc_ratios, labels, cols2 = [], [], []
for cohort, mut, wt, lab in idh_pairs:
    if len(mut) < 5 or len(wt) < 5:
        continue
    ratio = np.median(wt) / np.median(mut) if np.median(mut) > 0 else np.nan
    fc_ratios.append(ratio)
    labels.append(lab)
    cols2.append("#A32D2D" if cohort == "CGGA" else "#D9534F")
ax2.bar(range(len(fc_ratios)), fc_ratios, color=cols2, alpha=0.85,
        edgecolor="white")
ax2.axhline(1.0, color="#666", lw=0.8, ls="--")
ax2.set_xticks(range(len(labels)))
ax2.set_xticklabels(labels, fontsize=10)
ax2.set_ylabel("median ZP3\nIDH-WT / IDH-Mutan"
               "t")
ax2.set_title("IDH-WT vs Mutant\nZP3 ratio", fontsize=11)
for i, v in enumerate(fc_ratios):
    ax2.text(i, v * 1.03, f"{v:.2f}x", ha="center", va="bottom", fontsize=10,
             fontweight="bold")

plt.tight_layout()
fig.savefig(os.path.join(OUT, "tcga_subtype_direction_chain.png"),
            dpi=300, bbox_inches="tight")

# ---------- 6. 控制台结论 ----------
print("=== 亚型方向链值 ===")
print(val_df.to_string(index=False))
print()
print("=== IDH-WT vs Mutant 倍数 ===")
for cohort, mut, wt, lab in idh_pairs:
    if len(mut) < 5 and len(wt) < 5:
        continue
    if len(mut) >= 5 and len(wt) >= 5:
        r = np.median(wt) / np.median(mut)
        print(f"  {lab}: mut_med={np.median(mut):.3f}, wt_med={np.median(wt):.3f}, "
              f"ratio={r:.2f}x, (mut_n={len(mut)}, wt_n={len(wt)})")
print("\n已保存: tcga_subtype_direction_chain.png + subtype_direction_chain_values.csv")
