#!/usr/bin/env python
"""Corrected ZP3-immune marker correlation baseline for Article 1 (TCGA GBM/LGG).

AUDIT FIX (2026-08-24, mirrors the A2 GSE91061 Entrez correction):
    zp3_isoform_immune_real.py used ENTREZ_ZP3 = 8277, which maps to SP5, NOT
    ZP3 (correct Entrez = 7784). The frozen table zp3_immune_correlation_real.csv
    therefore quantified SP5-ZP3 correlations. Verified: LGG TREM2 rho=0.399,
    p=1.137e-21 reproduces exactly with 8277; with correct 7784, LGG
    TREM2 rho=0.164.

This script recomputes all 30 pre-specified immune/glial markers against
correct ZP3 (Entrez 7784) on the same TCGA GBM (n=166) / LGG (n=530)
cBioPortal RNA-seq V2 RSEM cohorts, applies BH-FDR within each study, and
writes a corrected baseline frozen table:
    article1/results/a1_zp3_immune_correlation_entrez7784.csv

It also writes a side-by-side comparison against the legacy (8277-based)
table for transparency:
    article1/results/a1_entrez_correction_comparison.csv
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

import time
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import requests

CBIO_API_BASE = "https://www.cbioportal.org/api"

ENTREZ_ZP3_CORRECT = 7784   # ZP3 (audited; 8277 = SP5, legacy error)

IMMUNE_GENES = {
    "TREM2": 54209, "CD68": 968, "CD163": 9332, "MRC1": 4360, "CD14": 929,
    "LYZ": 4069, "CSF1R": 1436, "ITGAM": 3684, "IL10": 3586, "TNF": 7124,
    "GZMA": 3001, "PRF1": 5551, "IFNG": 3458, "PDCD1": 5133, "CD274": 29126,
    "CTLA4": 1493, "LAG3": 3902, "FOXP3": 50943, "GPX4": 2879, "VEGFA": 7422,
    "CD8A": 925, "CD4": 920, "PTPRC": 5788, "AIF1": 199, "C1QA": 712,
    "TYROBP": 7305, "APOE": 348, "GFAP": 2670, "OLIG2": 10215, "SOX10": 6663,
}

STUDIES = {
    "GBM": {"profile": "gbm_tcga_rna_seq_v2_mrna", "sample_list": "gbm_tcga_rna_seq_v2_mrna"},
    "LGG": {"profile": "lgg_tcga_rna_seq_v2_mrna", "sample_list": "lgg_tcga_rna_seq_v2_mrna"},
}

LEGACY_CSV = os.path.join(RESULTS_DIR, "zp3_immune_correlation_real.csv")
OUT_CORRECTED = os.path.join(RESULTS_DIR, "a1_zp3_immune_correlation_entrez7784.csv")
OUT_COMPARE = os.path.join(RESULTS_DIR, "a1_entrez_correction_comparison.csv")


def cbio_get(endpoint, params=None):
    url = f"{CBIO_API_BASE}/{endpoint}"
    r = requests.get(url, params=params, timeout=60, headers={"Accept": "application/json"})
    r.raise_for_status()
    return r.json()


def fetch_gene(profile_id, sample_list, eid):
    data = cbio_get(
        f"molecular-profiles/{profile_id}/molecular-data",
        {"entrezGeneId": eid, "sampleListId": sample_list},
    )
    s = pd.Series({d["sampleId"]: d["value"] for d in data}, dtype=float)
    time.sleep(0.25)
    return s


def main():
    print("=" * 64)
    print("Corrected ZP3 (Entrez 7784) x 30-marker correlation, TCGA GBM/LGG")
    print("=" * 64)

    rows = []
    expr_cache = {}
    for study, info in STUDIES.items():
        zp3 = fetch_gene(info["profile"], info["sample_list"], ENTREZ_ZP3_CORRECT)
        print(f"[{study}] ZP3(7784): n={len(zp3)}")
        genes = {}
        for gname, eid in IMMUNE_GENES.items():
            genes[gname] = fetch_gene(info["profile"], info["sample_list"], eid)
        expr_cache[study] = pd.DataFrame(genes)

        for gname, gs in genes.items():
            mask = zp3.notna() & gs.notna()
            nz = (gs[mask] > 0).sum()
            if mask.sum() < 10 or nz < 10:
                continue
            rho, p = stats.spearmanr(zp3[mask], gs[mask])
            rows.append({"study": study, "gene": gname, "spearman_rho": rho,
                         "p_value": p, "n": int(mask.sum())})

    corr = pd.DataFrame(rows)
    corr["FDR"] = corr.groupby("study")["p_value"].transform(
        lambda s: multipletests(s, method="fdr_bh")[1])
    corr = corr.sort_values(["study", "spearman_rho"], ascending=[True, False])
    corr.to_csv(OUT_CORRECTED, index=False)
    print(f"\nSaved corrected baseline: {OUT_CORRECTED} ({len(corr)} rows)")

    # Side-by-side comparison vs legacy 8277-based table
    try:
        legacy = pd.read_csv(LEGACY_CSV)
        m = legacy.rename(columns={"spearman_rho": "rho_legacy_8277",
                                   "p_value": "p_legacy_8277",
                                   "FDR": "FDR_legacy_8277"})
        comp = corr.merge(m[["study", "gene", "rho_legacy_8277", "p_legacy_8277",
                             "FDR_legacy_8277"]], on=["study", "gene"], how="outer")
        comp["delta_rho"] = comp["spearman_rho"] - comp["rho_legacy_8277"]
        comp.to_csv(OUT_COMPARE, index=False)
        print(f"Saved comparison: {OUT_COMPARE}")

        print("\n=== Largest changes (|delta_rho| top 10) ===")
        top = comp.reindex(comp.delta_rho.abs().sort_values(ascending=False).index).head(10)
        for _, r in top.iterrows():
            print(f"  {r.study} {r.gene}: {r.rho_legacy_8277:.3f} -> {r.spearman_rho:.3f} "
                  f"(delta {r.delta_rho:+.3f})")

        print("\n=== Corrected headline markers ===")
        for study in ["GBM", "LGG"]:
            sub = corr[(corr.study == study) &
                       corr.gene.isin(["TREM2", "CSF1R", "C1QA", "CD274", "GPX4"])]
            for _, r in sub.iterrows():
                print(f"  {study} {r.gene}: rho={r.spearman_rho:.4f} "
                      f"p={r.p_value:.3e} FDR={r.FDR:.3e}")
        pos_lgg = corr[(corr.study == "LGG")]
        print(f"\nLGG positive correlations: {(pos_lgg.spearman_rho > 0).sum()}/"
              f"{len(pos_lgg)}; FDR<0.05: {(pos_lgg.FDR < 0.05).sum()}/{len(pos_lgg)}")
    except FileNotFoundError:
        print("Legacy table missing; comparison skipped.")

    print("\nNOTE: legacy zp3_immune_correlation_real.csv is RETAINED for audit trail "
          "but is INVALID as a ZP3 result (quantified SP5). Downstream analyses must "
          "use a1_zp3_immune_correlation_entrez7784.csv.")


if __name__ == "__main__":
    main()
