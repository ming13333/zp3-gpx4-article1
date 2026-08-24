#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
#38 Improvement 3: TCGA molecular subtype stratification—verify that ZP3 is highest in Mesenchymal (strongest immune infiltration)
=====================================================================
Data:
  · GBM subtypes (Verhaak 2010): gbm_tcga_pub2013 EXPRESSION_SUBTYPE
    (Mesenchymal / Classical / Proneural / Neural / G-CIMP) or
    lgggbm_tcga_pub TRANSCRIPTOME_SUBTYPE (ME / CL / PN / NE)
  · LGG subtypes (Cell 2016): lgggbm_tcga_pub IDH_1P19Q_SUBTYPE (Codel / Non-codel)
    + IDH_STATUS (WT / Mutant)
  · Local expression matrices: HiSeq_TCGA_gene.xena.gz from h2_bulk/TCGA.GBM.sampleMap and TCGA.LGG.sampleMap
    (log2 TPM, 23 genes x 153/509 samples)

Analysis:
  1. GBM: group ZP3 expression by transcriptome subtype (ME vs others) → verify ME is highest
  2. GBM: ZP3-immune signature association stratified by subtype
  3. LGG: group ZP3 expression by IDH status and 1p/19q
  4. LGG: ZP3-immune signature association stratified by IDH-codel

Products: tcga_subtype_zp3_*.csv + fig_tcga_subtype_zp3.png
"""
import os, sys, json, gzip, time
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
H2_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "output", "h2_bulk")  # output/h2_bulk
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article1", "results", "tcga_subtype_zp3_results")
os.makedirs(OUT, exist_ok=True)

M2_GENES = ["MRC1", "CD163", "MSR1", "ARG1", "TGFB1", "IL10", "VSIG4"]
TREG_GENES = ["FOXP3", "IL2RA", "CTLA4", "TIGIT"]
CHECKPT_GENES = ["CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3"]

CBIO = "https://www.cbioportal.org/api"


def fetch_clinical(study, attr, data_type="SAMPLE"):
    """Fetch cBioPortal clinical data, returns {sampleId: value}."""
    url = (f"{CBIO}/studies/{study}/clinical-data"
           f"?clinicalDataType={data_type}&attributeId={attr}&pageSize=10000")
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        records = json.loads(r.read().decode())
    return {rec["sampleId"]: rec.get("value", None) for rec in records}


def load_expression(cancer):
    """Load local Xena-format expression matrix. Returns DataFrame: index=sample, columns=gene."""
    p = os.path.join(H2_DIR, f"TCGA.{cancer}.sampleMap", "HiSeq_TCGA_gene.xena.gz")
    df = pd.read_csv(p, sep="\t", index_col=0, compression="gzip").T
    df.index.name = "sample"
    return df


def spearman_p(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return np.nan, np.nan, int(m.sum())
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), int(m.sum())


def score_zs(genes, df):
    """z-score consensus: standardize each gene across samples, then take the mean. df rows are samples, columns are genes."""
    avail = [g for g in genes if g in df.columns]
    if not avail:
        return pd.Series(np.nan, index=df.index)
    sub = df[avail].astype(float).T       # genes x samples
    v = sub.std(axis=1) > 0
    if not v.any():
        return pd.Series(np.nan, index=df.index)
    z = ((sub.loc[v] - sub.loc[v].mean(axis=1).values[:, None])
         / sub.loc[v].std(axis=1).values[:, None])
    return z.mean(axis=0)


def plain_mean(genes, df):
    """Plain mean score of a gene set, used for immune scores (log2 TPM scale)."""
    avail = [g for g in genes if g in df.columns]
    if not avail:
        return pd.Series(np.nan, index=df.index)
    return df[avail].astype(float).mean(axis=1)


def main():
    print("=== #38 TCGA molecular subtype stratification — ZP3/immune association ===\n")

    # ---- Step 1: Fetch subtypes ----
    print("1. Fetching cBioPortal molecular subtypes...")
    gbm_subtype = fetch_clinical("gbm_tcga_pub2013", "EXPRESSION_SUBTYPE")
    gbm_idh     = fetch_clinical("gbm_tcga_pub2013", "IDH1_MUTATION")
    lgg_idh1p19q = fetch_clinical("lgggbm_tcga_pub", "IDH_1P19Q_SUBTYPE")
    lgg_idh_stat = fetch_clinical("lgggbm_tcga_pub", "IDH_STATUS")
    lgg_tx_subtype = fetch_clinical("lgggbm_tcga_pub", "TRANSCRIPTOME_SUBTYPE")
    print(f"   GBM subtypes: {len(gbm_subtype)} | LGG IDH-1p19q: {len(lgg_idh1p19q)} | "
          f"LGG IDH: {len(lgg_idh_stat)} | LGG TX: {len(lgg_tx_subtype)}")

    # ---- Step 2: Load expression ----
    print("2. Loading local expression matrices...")
    expr_gbm = load_expression("GBM")
    expr_lgg = load_expression("LGG")
    print(f"   GBM: {expr_gbm.shape} | LGG: {expr_lgg.shape}")

    # ---- Step 3: GBM subtype analysis ----
    print("\n3. GBM molecular subtype analysis...")
    expr_gbm["subtype"] = expr_gbm.index.map(gbm_subtype.get)
    expr_gbm["idh"]     = expr_gbm.index.map(gbm_idh.get)
    gbm_typed = expr_gbm.dropna(subset=["subtype"])
    # Unified subtype grouping definition
    subtype_order = ["Mesenchymal", "Classical", "Proneural", "Neural"]
    gbm_typed["subtype_clean"] = gbm_typed["subtype"].apply(
        lambda x: x if x in subtype_order else "Other")
    print(f"   GBM samples with subtype: {len(gbm_typed)}")

    # A) ZP3 by subtype
    gbm_group_stats = []
    for st in subtype_order:
        grp = gbm_typed[gbm_typed["subtype"] == st]
        if len(grp) < 5:
            continue
        gbm_group_stats.append({
            "Cancer": "GBM", "Subtype": st,
            "N": len(grp),
            "ZP3_median": round(grp["ZP3"].median(), 4),
            "ZP3_mean": round(grp["ZP3"].mean(), 4)})

    # Mesenchymal vs others
    mes = gbm_typed[gbm_typed["subtype"] == "Mesenchymal"]["ZP3"].values
    non_mes = gbm_typed[gbm_typed["subtype"] != "Mesenchymal"]["ZP3"].values
    if len(mes) >= 5 and len(non_mes) >= 5:
        u, p_mes = stats.mannwhitneyu(mes, non_mes, alternative="two-sided")
        print(f"   Mesenchymal vs Others: median {np.median(mes):.3f} vs "
              f"{np.median(non_mes):.3f}, MWU p={p_mes:.4f}")

    # B) ZP3-immune association by subtype
    immuno_genes_all = sorted(set(M2_GENES + TREG_GENES + CHECKPT_GENES))
    rec_immuno = []
    for st in subtype_order:
        grp = gbm_typed[gbm_typed["subtype"] == st]
        if len(grp) < 15:
            continue
        for gset_name, genes in [("M2", M2_GENES), ("Treg", TREG_GENES),
                                  ("Checkpoint", CHECKPT_GENES)]:
            sc = plain_mean(genes, grp)
            rho, p, n = spearman_p(grp["ZP3"].values, sc.values)
            rec_immuno.append({
                "Cancer": "GBM", "Subtype": st, "GeneSet": gset_name,
                "Rho": rho, "P": p, "N": n, "Significant": bool(p < 0.05)})
    immuno_df = pd.DataFrame(rec_immuno)
    if not immuno_df.empty:
        print(f"   Immune associations (by subtype): {len(immuno_df)} records, "
              f"Sig: {immuno_df['Significant'].sum()}")

    # ---- Step 4: LGG subtype analysis ----
    print("\n4. LGG molecular subtype analysis...")
    expr_lgg["idh_1p19q"] = expr_lgg.index.map(lgg_idh1p19q.get)
    expr_lgg["idh_status"] = expr_lgg.index.map(lgg_idh_stat.get)
    expr_lgg["tx_subtype"] = expr_lgg.index.map(lgg_tx_subtype.get)
    lgg_typed = expr_lgg.dropna(subset=["idh_1p19q"])

    # IDH-codel vs Non-codel
    lgg_codel = lgg_typed[lgg_typed["idh_1p19q"] == "Codel"]
    lgg_noncodel = lgg_typed[lgg_typed["idh_1p19q"] == "Non-codel"]
    lgg_wt = lgg_typed[lgg_typed["idh_status"] == "WT"]
    lgg_mut = lgg_typed[lgg_typed["idh_status"] == "Mutant"]

    for label, grp in [("Codel", lgg_codel), ("Non-codel", lgg_noncodel),
                        ("IDH-WT", lgg_wt), ("IDH-Mut", lgg_mut)]:
        if len(grp) < 5:
            continue
        gbm_group_stats.append({
            "Cancer": "LGG", "Subtype": label,
            "N": len(grp),
            "ZP3_median": round(grp["ZP3"].median(), 4),
            "ZP3_mean": round(grp["ZP3"].mean(), 4)})

    # Codel vs Non-codel test
    if len(lgg_codel) >= 5 and len(lgg_noncodel) >= 5:
        u_c, p_c = stats.mannwhitneyu(lgg_codel["ZP3"].values,
                                       lgg_noncodel["ZP3"].values,
                                       alternative="two-sided")
        print(f"   LGG Codel vs Non-codel: median {lgg_codel['ZP3'].median():.3f} vs "
              f"{lgg_noncodel['ZP3'].median():.3f}, MWU p={p_c:.4f}")

    # LGG TX subtype ZP3
    lgg_tx = lgg_typed.dropna(subset=["tx_subtype"])
    for st in ["ME", "CL", "PN", "NE"]:
        grp = lgg_tx[lgg_tx["tx_subtype"] == st]
        if len(grp) < 5:
            continue
        gbm_group_stats.append({
            "Cancer": "LGG_TX", "Subtype": st,
            "N": len(grp),
            "ZP3_median": round(grp["ZP3"].median(), 4),
            "ZP3_mean": round(grp["ZP3"].mean(), 4)})

    # LGG immune association by codel
    for label, grp in [("Codel", lgg_codel), ("Non-codel", lgg_noncodel)]:
        if len(grp) < 15:
            continue
        for gset_name, genes in [("M2", M2_GENES), ("Treg", TREG_GENES),
                                  ("Checkpoint", CHECKPT_GENES)]:
            sc = plain_mean(genes, grp)
            rho, p, n = spearman_p(grp["ZP3"].values, sc.values)
            rec_immuno.append({
                "Cancer": "LGG", "Subtype": label, "GeneSet": gset_name,
                "Rho": rho, "P": p, "N": n, "Significant": bool(p < 0.05)})

    # ---- Step 5: Save CSV ----
    pd.DataFrame(gbm_group_stats).to_csv(
        os.path.join(OUT, "tcga_subtype_zp3_expression.csv"), index=False)
    pd.DataFrame(rec_immuno).to_csv(
        os.path.join(OUT, "tcga_subtype_zp3_immuno.csv"), index=False)
    gbm_typed[["ZP3", "subtype", "subtype_clean", "idh"]].to_csv(
        os.path.join(OUT, "gbm_zp3_subtype_annotated.csv"))
    lgg_typed[["ZP3", "idh_1p19q", "idh_status", "tx_subtype"]].to_csv(
        os.path.join(OUT, "lgg_zp3_subtype_annotated.csv"))

    # ---- Step 6: Figures ----
    print("\n5. Generating molecular subtype figures...")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

    # a) GBM subtype ZP3
    ax = axes[0]
    colors = {"Mesenchymal": "#A32D2D", "Classical": "#378ADD",
              "Proneural": "#1D9E75", "Neural": "#888780"}
    order = [s for s in subtype_order if s in gbm_typed["subtype"].values]
    bp_data = [gbm_typed[gbm_typed["subtype"] == s]["ZP3"].dropna().values
               for s in order]
    bp = ax.boxplot(bp_data, labels=order, patch_artist=True, showfliers=False)
    for patch, s in zip(bp["boxes"], order):
        patch.set_facecolor(colors.get(s, "#CCCCCC"))
    ax.set_ylabel("ZP3 log2(TPM)")
    ax.set_title(f"GBM ZP3 by Verhaak Subtype\n(n={len(gbm_typed)})", fontsize=11)

    # b) LGG IDH-codel ZP3
    ax = axes[1]
    codel_labels = ["Codel", "Non-codel"]
    codel_data = [lgg_codel["ZP3"].dropna().values,
                  lgg_noncodel["ZP3"].dropna().values]
    bp2 = ax.boxplot(codel_data, labels=codel_labels, patch_artist=True,
                     showfliers=False)
    bp2["boxes"][0].set_facecolor("#1D9E75")
    bp2["boxes"][1].set_facecolor("#888780")
    ax.set_ylabel("ZP3 log2(TPM)")
    ax.set_title(f"LGG ZP3 by 1p/19q\n(n={len(lgg_typed)})", fontsize=11)

    # c) Subtype ZP3 median bar chart
    ax = axes[2]
    stats_df = pd.DataFrame(gbm_group_stats)
    stats_df["label"] = stats_df["Cancer"] + ":" + stats_df["Subtype"]
    stats_df = stats_df.sort_values("ZP3_median", ascending=True)
    ax.barh(range(len(stats_df)), stats_df["ZP3_median"].values,
            color=[colors.get(s, "#888780")
                   for s in stats_df["Subtype"].values])
    ax.set_yticks(range(len(stats_df)))
    ax.set_yticklabels(stats_df["label"], fontsize=9)
    ax.set_xlabel("median ZP3 log2(TPM)")
    ax.set_title("ZP3 across molecular subtypes", fontsize=11)

    plt.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_tcga_subtype_zp3.png"),
                dpi=300, bbox_inches="tight")
    print("   Saved fig_tcga_subtype_zp3.png")

    # ---- Conclusion ----
    print("\n=== #38 Conclusion ===")
    print("  Testing hypothesis: Is ZP3 highest in Mesenchymal (strongest immune infiltration) GBM?")
    if "ZP3_median" in pd.DataFrame(gbm_group_stats).columns:
        mes_row = [r for r in gbm_group_stats
                   if r["Cancer"] == "GBM" and r["Subtype"] == "Mesenchymal"]
        if mes_row:
            print(f"  GBM Mesenchymal: ZP3 median = {mes_row[0]['ZP3_median']}, "
                  f"n = {mes_row[0]['N']}")
            others = [r for r in gbm_group_stats
                      if r["Cancer"] == "GBM" and r["Subtype"] != "Mesenchymal"]
            others_med = np.median([r["ZP3_median"] for r in others])
            print(f"  Others median = {others_med:.3f}, "
                  f"ratio = {mes_row[0]['ZP3_median']/others_med:.2f}x")
    print(f"  LGG Codel n={len(lgg_codel)} vs Non-codel n={len(lgg_noncodel)}")
    print("\n=== #38 Done ===")


if __name__ == "__main__":
    main()
