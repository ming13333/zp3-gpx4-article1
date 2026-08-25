#!/usr/bin/env python
"""Ferroptosis covariate sensitivity analysis for Article 1 (A1).

Question addressed (per codex-auto-review / gpt-5.4-mini strategic review, 2026-08-24):
    Is the ZP3-immune association in TCGA GBM/LGG explained by ferroptosis state?
Design:
    - Pull a ferroptosis gene panel (8 genes) from cBioPortal for the SAME
      TCGA GBM (n=166) / LGG (n=530) RNA-seq v2 RSEM cohorts used by
      zp3_isoform_immune_real.py, so sample sets are identical.
    - Build a per-sample ferroptosis score = mean of gene-wise z-scores of
      log2(RSEM + 1).
    - For each immune feature significant at baseline (FDR < 0.05 in the frozen
      zp3_immune_correlation_real.csv), compute partial Spearman correlation
      rho(ZP3, immune | ferroptosis score) via rank-transform residualization.
    - Report rho_base -> rho_adj, delta-rho, partial-correlation p-value and a
      Fisher-z 95% CI for rho_adj.
Boundary statement (mandatory):
    This is an ALTERNATIVE-EXPLANATION sensitivity analysis. It tests whether
    the association survives adjustment; it does NOT establish that ZP3 acts
    independently of ferroptosis, nor any causal claim.

Output: article1/results/a1_ferroptosis_sensitivity.csv (frozen table)
"""
import os as _os

def _project_root():
    d = _os.path.dirname(_os.path.abspath(__file__))
    while True:
        if _os.path.isdir(_os.path.join(d, "output")):
            return d
        p = _os.path.dirname(d)
        if p == d:
            break
        d = p
    return _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

ROOT = _project_root()
import os
RESULTS_DIR = os.path.join(ROOT, "article1", "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

import requests
import time
import numpy as np
import pandas as pd
from scipy import stats

CBIO_API_BASE = "https://www.cbioportal.org/api"

# Ferroptosis panel (codex-auto-review recommendation, 2026-08-24)
# Entrez IDs verified against NCBI Gene.
FERROPTOSIS_GENES = {
    "GPX4": 2879,
    "SLC7A11": 23657,
    "ACSL4": 2182,
    "AIFM2": 137275,   # aka FSP1
    "TFRC": 7037,
    "FTH1": 2495,
    "NCOA4": 8030,
    "ALOX15": 246,
}

STUDIES = {
    "GBM": {"profile": "gbm_tcga_rna_seq_v2_mrna", "sample_list": "gbm_tcga_rna_seq_v2_mrna"},
    "LGG": {"profile": "lgg_tcga_rna_seq_v2_mrna", "sample_list": "lgg_tcga_rna_seq_v2_mrna"},
}

BASELINE_CSV = os.path.join(RESULTS_DIR, "a1_zp3_immune_correlation_entrez7784.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "a1_ferroptosis_sensitivity.csv")

FDR_BASELINE = 0.05   # which immune features to re-test


def cbio_get(endpoint, params=None):
    url = f"{CBIO_API_BASE}/{endpoint}"
    r = requests.get(url, params=params, timeout=60, headers={"Accept": "application/json"})
    r.raise_for_status()
    return r.json()


def fetch_expression(profile_id, sample_list, entrez_map):
    """Fetch expression rows for all genes; return wide DataFrame."""
    frames = []
    for name, eid in entrez_map.items():
        data = cbio_get(
            f"molecular-profiles/{profile_id}/molecular-data",
            {"entrezGeneId": eid, "sampleListId": sample_list},
        )
        s = pd.Series({d["sampleId"]: d["value"] for d in data}, name=name)
        frames.append(s)
        print(f"    {name}: {len(s)} samples")
        time.sleep(0.3)
    df = pd.concat(frames, axis=1)
    return df


def fisher_ci(rho, n, alpha=0.05):
    """95% CI for a Pearson-style correlation via Fisher z (approximation)."""
    if not np.isfinite(rho) or abs(rho) >= 1 or n < 4:
        return (np.nan, np.nan)
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - 3)
    zc = stats.norm.ppf(1 - alpha / 2)
    return (np.tanh(z - zc * se), np.tanh(z + zc * se))


def partial_spearman(x, y, z):
    """Partial Spearman correlation of x,y controlling z (rank-residual method)."""
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    rz = pd.Series(z).rank().values

    def resid(a, b):
        b1 = np.column_stack([b, np.ones_like(b)])
        coef, *_ = np.linalg.lstsq(b1, a, rcond=None)
        return a - b1 @ coef

    ex = resid(rx, rz)
    ey = resid(ry, rz)
    if np.std(ex) == 0 or np.std(ey) == 0:
        return np.nan, np.nan
    r, p = stats.pearsonr(ex, ey)
    return r, p


def main():
    print("=" * 64)
    print("A1 ferroptosis-covariate sensitivity analysis (TCGA GBM/LGG)")
    print("=" * 64)

    # 1) Load frozen baseline correlations
    base = pd.read_csv(BASELINE_CSV)
    print(f"Baseline frozen table: {len(base)} rows "
          f"(GBM {sum(base.study=='GBM')}, LGG {sum(base.study=='LGG')})")

    # 2) Fetch ferroptosis panel per study
    expr = {}
    for study, info in STUDIES.items():
        print(f"\nFetching ferroptosis panel for {study} ...")
        df = fetch_expression(info["profile"], info["sample_list"], FERROPTOSIS_GENES)
        expr[study] = df
        print(f"  {study}: matrix {df.shape}")

    # Also need ZP3 + immune genes on the same samples.
    IMMUNE_ENTREZ = {
        "TREM2": 54209, "CD68": 968, "CD163": 9332, "MRC1": 4360, "CD14": 929,
        "LYZ": 4069, "CSF1R": 1436, "ITGAM": 3684, "IL10": 3586, "TNF": 7124,
        "GZMA": 3001, "PRF1": 5551, "IFNG": 3458, "PDCD1": 5133, "CD274": 29126,
        "CTLA4": 1493, "LAG3": 3902, "FOXP3": 50943, "VEGFA": 7422, "CD8A": 925,
        "CD4": 920, "PTPRC": 5788, "AIF1": 199, "C1QA": 712, "TYROBP": 7305,
        "APOE": 348, "GFAP": 2670, "OLIG2": 10215, "SOX10": 6663,
        "GPX4": 2879,  # in ferroptosis panel too; needed here as an outcome gene
    }
    zp3_immune = {}
    for study, info in STUDIES.items():
        print(f"\nFetching ZP3 + immune genes for {study} ...")
        full_map = dict(IMMUNE_ENTREZ); full_map["ZP3"] = 7784  # audited correct ZP3 (8277=SP5 legacy error)
        # NOTE: baseline table a1_zp3_immune_correlation_entrez7784.csv was
        # recomputed with Entrez 7784; this script uses the same ID.
        df = fetch_expression(info["profile"], info["sample_list"], full_map)
        zp3_immune[study] = df
        print(f"  {study}: matrix {df.shape}")

    results = []
    for study in STUDIES:
        fer = expr[study]
        zi = zp3_immune[study]

        # Align to intersection of samples (should be ~identical)
        common = fer.index.intersection(zi.index)
        fer = fer.loc[common]; zi = zi.loc[common]
        n = len(common)

        # Ferroptosis score: mean of per-gene z-scores of log2(x+1)
        lg = np.log2(fer.astype(float).clip(lower=0) + 1)
        zmat = lg.apply(lambda col: (col - col.mean()) / col.std(ddof=0))
        fer_score = zmat.mean(axis=1)

        # Baseline sanity: reproduce rho(ZP3, GPX4) from frozen table
        zp3v = zi["ZP3"].astype(float)
        gpx4_check = stats.spearmanr(zp3v, fer["GPX4"].astype(float))

        sig_features = base[(base.study == study) & (base.FDR < FDR_BASELINE)]["gene"].unique()

        print(f"\n[{study}] n={n}; ferroptosis score mean={fer_score.mean():.3f} sd={fer_score.std():.3f}")
        print(f"  GPX4-ZP3 baseline check: rho={gpx4_check.statistic:.4f} p={gpx4_check.pvalue:.3g}")
        print(f"  Immune features at baseline FDR<{FDR_BASELINE}: {len(sig_features)}")

        for gene in sorted(sig_features):
            yv = zi[gene].astype(float)
            mask = zp3v.notna() & yv.notna() & fer_score.notna()
            x_, y_, z_ = zp3v[mask], yv[mask], fer_score[mask]

            rb = stats.spearmanr(x_, y_)
            rp, pp = partial_spearman(x_, y_, z_)
            ci_lo, ci_hi = fisher_ci(rp, int(mask.sum()))

            results.append({
                "study": study,
                "gene": gene,
                "n": int(mask.sum()),
                "rho_base": round(rb.statistic, 4),
                "rho_adj": None if not np.isfinite(rp) else round(rp, 4),
                "delta_rho": None if not np.isfinite(rp) else round(rp - rb.statistic, 4),
                "ci95_low": None if not np.isfinite(ci_lo) else round(ci_lo, 4),
                "ci95_high": None if not np.isfinite(ci_hi) else round(ci_hi, 4),
                "p_partial": None if not np.isfinite(pp) else f"{pp:.3e}",
                "direction_retained": None if not np.isfinite(rp) else bool(np.sign(rb.statistic) == np.sign(rp)),
            })

    out = pd.DataFrame(results)
    out.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved frozen table: {OUTPUT_CSV} ({len(out)} rows)")

    # Summary statistics
    print("\n=== SUMMARY ===")
    for study in STUDIES:
        sub = out[out.study == study].dropna(subset=["rho_adj"])
        ret = sub.direction_retained.mean() * 100
        med_delta = sub.delta_rho.median()
        still_sig = sum(
            float(r.p_partial) < 0.05 for _, r in sub.iterrows()
        )
        print(f"{study}: {len(sub)} features adjusted; "
              f"direction retained {ret:.0f}%; median delta_rho={med_delta:+.4f}; "
              f"partial p<0.05 in {still_sig}/{len(sub)}")

    print("\nBoundary note: this analysis tests whether associations survive "
          "adjustment for a ferroptosis-state score; it does NOT establish that "
          "ZP3 is independent of ferroptosis biology or imply causality.")


if __name__ == "__main__":
    main()
