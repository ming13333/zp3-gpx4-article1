# -*- coding: utf-8 -*-
"""
H2 (prognostic association) + H3 (immunosuppressive TME association) — TCGA GBM/LGG bulk RNA-seq
Data: HiSeq_TCGA_gene.xena.gz downloaded from UCSC Xena (expression, rows=genes, columns=samples)
      GBM/LGG_clinicalMatrix.gz (clinical, containing OS_MONTHS / OS_STATUS)
Descriptive-level analysis, honestly labeled.
"""
import os, gzip, sys, numpy as np, pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))

# Use the log-rank validated against lifelines/standard implementations (eliminates two inconsistent legacy hand-written implementations)
sys.path.insert(0, os.path.join(BASE, "..", "common"))
from stats_utils import logrank  # noqa: E402

# ---- Immunosuppression-related marker genes (common in literature, for descriptive association) ----
IMMUNOSUPP_GENES = ["TGFB1", "IL10", "FOXP3", "CD274", "PDCD1", "CTLA4",
                    "MRC1", "CD163", "VSIG4", "ARG1", "IDO1", "VEGFA",
                    "CCL2", "CXCL12", "MSR1", "TREM2"]
# M2 / TAM biased
M2_GENES = ["MRC1", "CD163", "MSR1", "ARG1", "TGFB1", "IL10", "VSIG4"]
# Treg
TREG_GENES = ["FOXP3", "IL2RA", "CTLA4", "TIGIT"]
# Checkpoints
CHECKPT_GENES = ["CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3"]

def load_expr(path):
    """Xena expression matrix: rows=gene symbols, columns=samples. Returns DataFrame."""
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    return df

def load_clin(path):
    """Xena clinical matrix: try parsing as attributes x samples or samples x attributes."""
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    # If OS_MONTHS is in index -> attributes x samples
    if "OS_MONTHS" in df.index or "OS_MONTHS" in df.columns:
        if "OS_MONTHS" in df.index:
            return df.T  # -> samples x attributes
        return df
    return df

def extract_os(clin, sample_ids):
    """Extract OS_MONTHS / OS_STATUS from the clinical table, matching expression samples."""
    # clin: samples x attributes
    common = [s for s in sample_ids if s in clin.index]
    sub = clin.loc[common]
    # Find OS columns (flexible matching)
    os_time_col = next((c for c in sub.columns if c.upper() in ("OS_MONTHS", "_OS_MONTHS", "OS_MONTHS")), None)
    os_stat_col = next((c for c in sub.columns if "OS_STATUS" in c.upper() or c.upper() in ("_OS_STATUS",)), None)
    if os_time_col is None or os_stat_col is None:
        # Print available columns to help debugging
        print("  Available clinical columns (first 30):", list(sub.columns[:30]))
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
    print(f"Cohort: {name}")
    expr = load_expr(expr_path)
    clin = load_clin(clin_path)
    print(f"  Expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")
    if "ZP3" not in expr.index:
        print("  !! No ZP3 row in expression matrix, skipping")
        return
    zp3 = expr.loc["ZP3"]
    # Keep only numeric samples
    zp3 = pd.to_numeric(zp3, errors="coerce").dropna()
    osd = extract_os(clin, zp3.index.tolist())
    if osd is None or len(osd) < 20:
        print("  !! Insufficient OS data, skipping survival")
        return
    # Merge
    merged = pd.concat([zp3.rename("ZP3"), osd], axis=1).dropna()
    merged = merged[merged["time"] > 0]
    print(f"  Samples available for survival: {len(merged)}")
    # Median dichotomization
    med = merged["ZP3"].median()
    merged["group"] = (merged["ZP3"] > med).astype(int)
    hi = merged[merged.group == 1]; lo = merged[merged.group == 0]
    chi2, p = logrank(merged["time"].values, merged["event"].values, merged["group"].values)
    print(f"  ZP3 median={med:.3f} | High n={len(hi)} vs Low n={len(lo)}")
    print(f"  H2 logrank: chi2={chi2:.3f}, p={p:.4g}")
    # Direction: event rate in High group
    rate_hi = hi["event"].mean(); rate_lo = lo["event"].mean()
    print(f"  Event rate High={rate_hi:.2f} Low={rate_lo:.2f} -> {'High ZP3 worse prognosis' if rate_hi>rate_lo else 'High ZP3 better prognosis'}")

    # H3: association between ZP3 and immunosuppressive markers
    # Note: expression is right-skewed and non-normal, so Spearman rank correlation is used (originally Pearson, corrected)
    print(f"  --- H3: ZP3 vs immunosuppressive markers (Spearman rho) ---")
    rows = []
    h3_genes = list(dict.fromkeys(IMMUNOSUPP_GENES + M2_GENES + TREG_GENES + CHECKPT_GENES))  # deduplicate (same gene in multiple sets)
    for gene in h3_genes:
        if gene in expr.index:
            g = pd.to_numeric(expr.loc[gene], errors="coerce")
            gg = pd.concat([merged["ZP3"], g.rename(gene)], axis=1).dropna()
            if len(gg) > 20:
                r, pp = stats.spearmanr(gg["ZP3"], gg[gene])
                rows.append((gene, round(r, 3), round(pp, 4), len(gg)))
    h3 = pd.DataFrame(rows, columns=["gene", "spearman_rho", "p", "n"]).sort_values("spearman_rho", ascending=False)
    print(h3.to_string(index=False))
    # Composite immunosuppression score (z-mean of available M2+TREG+CHECKPT genes)
    sig = [g for g in M2_GENES + TREG_GENES + CHECKPT_GENES if g in expr.index]
    if sig:
        sub = expr.loc[sig, merged.index].apply(pd.to_numeric, errors="coerce")
        z = (sub - sub.mean()) / sub.std()
        immuno_score = z.mean(axis=0)
        gg = pd.concat([merged["ZP3"], immuno_score.rename("immuno_score")], axis=1).dropna()
        r, pp = stats.spearmanr(gg["ZP3"], gg["immuno_score"])
        print(f"  Composite immunosuppression score vs ZP3: rho={r:.3f}, p={pp:.4g}, n={len(gg)}")
    return merged, h3

if __name__ == "__main__":
    results = {}
    results["GBM"] = analyze_cohort("TCGA GBM",
        os.path.join(BASE, "TCGA.GBM.sampleMap/HiSeq_TCGA_gene.xena.gz"),
        os.path.join(BASE, "TCGA.GBM.sampleMap/GBM_clinicalMatrix.gz"))
    results["LGG"] = analyze_cohort("TCGA LGG",
        os.path.join(BASE, "TCGA.LGG.sampleMap/HiSeq_TCGA_gene.xena.gz"),
        os.path.join(BASE, "TCGA.LGG.sampleMap/LGG_clinicalMatrix.gz"))
    print("\n=== H2/H3 analysis complete (descriptive level, requires external validation) ===")
