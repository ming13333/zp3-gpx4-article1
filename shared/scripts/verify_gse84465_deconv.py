# -*- coding: utf-8 -*-
"""
Method 1 · Single-cell origin identification — GSE84465 independent cohort cell-level replication (cross-cohort robustness test)

Background: GSE182109 subset h5ad (10x 3', HVG space) in verify_sources_deconv.py Test A
  showed that 89.3% of ZP3+ myeloid cells co-express TREM2+ (OR=20.5, p=8.5e-11).
This script in the second independent GBM single-cell cohort GSE84465 (Darmanis 2017, SMART-seq2 full-length,
3589 cells, genes×cells, raw count scale) redoes the same test to evaluate cross-cohort robustness.

Outputs:
  - Co-enrichment table: all cells / within myeloid / MG-TAM-DC subsets (Fisher exact OR+p)
  - Logistic regression: ZP3 expression -> TREM2+ (per-unit OR + likelihood ratio p)
  - ZP3 source lineage analysis: are ZP3+ cells enriched in myeloid? Lineage marker characteristics of ZP3-high cells
  - Item-by-item comparison table with GSE182109 main cohort

2026-08-10  Craft mode execution (supplement cross-cohort robustness)
"""
import os, gzip
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.optimize import minimize

OUT = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(os.path.dirname(OUT), "h1_pilot", "GSE84465_GBM_All_data.csv.gz")
LOG = []

def log(msg=""):
    print(msg)
    LOG.append(str(msg))

# ---- marker gating (identical to h1_replicate_gse84465.py, ensuring comparability) ----
pan_myeloid = ["CD68", "LYZ", "C1QA", "C1QB", "ITGAM", "CSF1R", "CD14"]
MG  = ["CX3CR1", "P2RY12", "TMEM119", "SALL1", "SIGLEC11"]
TAM = ["CD163", "VSIG4", "MRC1", "MSR1", "FOLR2"]
DC  = ["CLEC9A", "FCER1A", "CD1C", "LAMP3", "BATF3", "ITGAX"]

def fisher_table(a, b):
    a = np.asarray(a, bool); b = np.asarray(b, bool)
    aa = int(((a) & (b)).sum()); ab = int(((a) & (~b)).sum())
    ba = int(((~a) & (b)).sum()); bb = int(((~a) & (~b)).sum())
    OR, p = stats.fisher_exact([[aa, ab], [ba, bb]])
    return aa, ab, ba, bb, OR, p

def logistic_or(x, y):
    """logit(y)=b0+b1*x -> (OR=exp(b1), b0, b1)"""
    x = np.asarray(x, float); y = np.asarray(y, float)
    Xd = np.column_stack([np.ones_like(x), x])
    def negll(b):
        z = Xd @ b
        p = 1/(1+np.exp(-z))
        p = np.clip(p, 1e-12, 1-1e-12)
        return -(y*np.log(p) + (1-y)*np.log(1-p)).sum()
    res = minimize(negll, np.zeros(2), method='BFGS')
    return float(np.exp(res.x[1])), float(res.x[0]), float(res.x[1])

def lrt_p(x, y, b0, b1):
    """Likelihood ratio test p (H0: b1=0)"""
    x = np.asarray(x, float); y = np.asarray(y, float)
    Xd = np.column_stack([np.ones_like(x), x])
    def ll(b):
        z = Xd @ b
        p = 1/(1+np.exp(-z))
        p = np.clip(p, 1e-12, 1-1e-12)
        return (y*np.log(p) + (1-y)*np.log(1-p)).sum()
    ll1 = ll([b0, b1])
    ll0 = ll([np.log(y.mean()/(1-y.mean())), 0.0])
    return float(stats.chi2.sf(2*(ll1 - ll0), 1))

def main():
    with gzip.open(PATH, "rt") as f:
        first = f.readline()
    sep = r"\s+" if first.count(" ") > first.count(",") else ","
    df = pd.read_csv(PATH, sep=sep, index_col=0, compression="gzip")
    expr = df.T.apply(pd.to_numeric, errors="coerce").fillna(0)

    log("=" * 74)
    log("GSE84465 cell-level ZP3↔TREM2 co-enrichment + ZP3 source lineage (independent cohort replication)")
    log("  cells=%d genes=%d | platform: SMART-seq2 full-length | index=plate position (no official annotation)"
        % expr.shape)
    zp3 = expr["ZP3"].values.astype(float)
    trem2 = expr["TREM2"].values.astype(float)
    zp3_pos = zp3 > 0; trem2_pos = trem2 > 0
    n = len(zp3)
    log("  background: ZP3+ =%d (%.2f%%) | TREM2+ =%d (%.2f%%)"
        % (zp3_pos.sum(), 100*zp3_pos.mean(), trem2_pos.sum(), 100*trem2_pos.mean()))
    log("  data scale check: GAPDH p50=%.0f (raw count magnitude) | ZP3 max=%.0f | TREM2 max=%.0f"
        % (np.median(expr["GAPDH"].values.astype(float)), zp3.max(), trem2.max()))

    # ---- marker gating ----
    def cmean(gl):
        ix = [g for g in gl if g in expr.columns]
        return expr[ix].mean(axis=1) if ix else pd.Series(0.0, index=expr.index)
    pm = cmean(pan_myeloid).values.astype(float)
    mg = cmean(MG); tam = cmean(TAM); dc = cmean(DC)
    myeloid = pm > 0
    sub = pd.DataFrame({"MG": mg, "TAM": tam, "DC": dc}).loc[myeloid]
    best = sub.idxmax(axis=1); val = sub.max(axis=1)
    subclass_pool = best.where(val > 0, "Unassigned")
    subclass_full = pd.Series("Not_myeloid", index=expr.index)
    subclass_full.loc[myeloid] = subclass_pool
    sc = subclass_full.values.astype(str)
    log("\ngating (consistent with GSE84465 H1 replication): pan-myeloid n=%d (%.1f%%) | TAM=%d MG=%d DC=%d Un=%d"
        % (myeloid.sum(), 100*myeloid.mean(),
           int((sc == "TAM").sum()), int((sc == "MG").sum()),
           int((sc == "DC").sum()), int((sc == "Unassigned").sum())))

    # ---- Table 1: co-enrichment ----
    log("\n[Table 1] Is ZP3+ enriched in TREM2+")
    rows = []
    for label, mask in [("all cells", np.ones(n, bool)), ("within myeloid", myeloid)]:
        aa, ab, ba, bb, OR, p = fisher_table(zp3_pos[mask], trem2_pos[mask])
        frac = aa/(aa+ab) if aa+ab else float('nan')
        bg = (trem2_pos[mask]).mean()
        log("  [%s] n=%d | ZP3+&TREM2+=%d  ZP3+&TREM2-=%d | TREM2+ proportion among ZP3+ =%.1f%% (background %.1f%%) | OR=%.2f p=%.3g"
            % (label, int(mask.sum()), aa, ab, 100*frac, 100*bg, OR, p))
        rows.append({"level": label, "n": int(mask.sum()), "n_zp3pos": aa+ab,
                     "n_zp3pos_trem2pos": aa, "frac_zp3pos_trem2pos": round(frac, 4),
                     "bg_trem2pos": round(float(bg), 4), "OR": round(OR, 2), "p": p})
    for s in ["MG", "TAM", "DC"]:
        mask = sc == s
        if int(mask.sum()) < 10:
            continue
        aa, ab, ba, bb, OR, p = fisher_table(zp3_pos[mask], trem2_pos[mask])
        frac = aa/(aa+ab) if aa+ab else float('nan')
        bg = (trem2_pos[mask]).mean()
        log("  [%s] n=%d | ZP3+&TREM2+=%d/%d | within-subset TREM2+ background=%.1f%% | OR=%.2f p=%.3g"
            % (s, int(mask.sum()), aa, aa+ab, 100*bg, OR, p))
        rows.append({"level": s, "n": int(mask.sum()), "n_zp3pos": aa+ab,
                     "n_zp3pos_trem2pos": aa,
                     "frac_zp3pos_trem2pos": (round(frac, 4) if frac == frac else None),
                     "bg_trem2pos": round(float(bg), 4), "OR": round(OR, 2), "p": p})

    # ---- Table 2: logistic regression ----
    log("\n[Table 2] ZP3 expression -> TREM2+ (log1p, likelihood ratio test)")
    lr_rows = []
    x_all = np.log1p(zp3); y_all = trem2_pos.astype(float)
    or_a, b0, b1 = logistic_or(x_all, y_all)
    p_a = lrt_p(x_all, y_all, b0, b1)
    log("  All cells: OR=%.2f / unit log1p ZP3, p=%.3g" % (or_a, p_a))
    lr_rows.append({"level": "All cells", "per_unit_OR": round(or_a, 2), "p": p_a})
    if myeloid.sum() > 0:
        x_m = np.log1p(zp3[myeloid]); y_m = trem2_pos[myeloid].astype(float)
        or_m, b0m, b1m = logistic_or(x_m, y_m)
        p_m = lrt_p(x_m, y_m, b0m, b1m)
        log("  Within myeloid: OR=%.2f / unit log1p ZP3, p=%.3g" % (or_m, p_m))
        lr_rows.append({"level": "Within myeloid", "per_unit_OR": round(or_m, 2), "p": p_m})

    # ---- Table 3: Whether ZP3 is myeloid-derived (direct test independent of TREM2)----
    log("\n[Table 3] Are ZP3+ cells enriched in myeloid (direct test of 'bulk ZP3 from myeloid?')")
    src_rows = []
    for pth in [0, 5]:
        a = zp3 > pth
        if int(a.sum()) < 10:
            continue
        aa, ab, ba, bb, OR, p = fisher_table(a, myeloid)
        frac = aa/(aa+ab) if aa+ab else float('nan')
        bg = myeloid.mean()
        log("  ZP3>%d: myeloid in ZP3+=%.1f%% (background %.1f%%) | OR=%.2f p=%.3g (n=%d)"
            % (pth, 100*frac, 100*bg, OR, p, int(a.sum())))
        src_rows.append({"zp3_threshold": pth, "n_zp3pos": int(a.sum()),
                         "frac_zp3pos_myeloid": round(frac, 4),
                         "bg_myeloid": round(float(bg), 4),
                         "OR": round(OR, 2), "p": p})

    # ---- Table 4: Lineage marker features of ZP3-high cells ----
    log("\n[Table 4] ZP3-high cells (ZP3>10, n=%d) vs ZP3=0 (n=%d) lineage marker log2FC"
        % (int((zp3 > 10).sum()), int((zp3 == 0).sum())))
    markers = {
        "Tumor/Neural": ["EGFR", "OLIG2", "SOX2", "GFAP", "NES", "TOP2A"],
        "Myeloid": ["CD68", "LYZ", "CD163", "CX3CR1", "P2RY12", "TMEM119"],
        "T/NK": ["NKG7", "CD3D"],
        "Endothelial/Pericyte": ["VWF", "PDGFRB", "COL1A1"],
    }
    hi_mask = (zp3 > 10); lo_mask = (zp3 == 0)
    marker_rows = []
    for grp, gl in markers.items():
        avail = [g for g in gl if g in expr.columns]
        if not avail:
            continue
        hi = expr.loc[hi_mask, avail].mean(axis=0)
        lo = expr.loc[lo_mask, avail].mean(axis=0)
        lfc = np.log2((hi.values + 1) / (lo.values + 1))
        for i, g in enumerate(avail):
            log("  %-8s (%-6s): hi=%.1f lo=%.1f log2FC=%+.2f"
                % (g, grp, hi.iloc[i], lo.iloc[i], lfc[i]))
            marker_rows.append({"gene": g, "group": grp,
                                "zp3hi_mean": round(float(hi.iloc[i]), 2),
                                "zp3lo_mean": round(float(lo.iloc[i]), 2),
                                "log2FC": round(float(lfc[i]), 2)})

    # ---- Save ----
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "gse84465_coenrichment.csv"), index=False)
    pd.DataFrame(lr_rows).to_csv(os.path.join(OUT, "gse84465_logistic.csv"), index=False)
    pd.DataFrame(src_rows).to_csv(os.path.join(OUT, "gse84465_source_myeloid.csv"), index=False)
    pd.DataFrame(marker_rows).to_csv(os.path.join(OUT, "gse84465_zp3hi_markers.csv"), index=False)
    log("\nSaved: gse84465_coenrichment.csv / gse84465_logistic.csv / gse84465_source_myeloid.csv / gse84465_zp3hi_markers.csv")
    log("=" * 74)

if __name__ == "__main__":
    main()
