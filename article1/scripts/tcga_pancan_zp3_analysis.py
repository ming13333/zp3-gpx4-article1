#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TCGA pan-cancer analysis of tissue specificity of ZP3 and immune features (real data version)
=================================================================
Data: UCSC Xena 'TCGA TARGET GTEx' RSEM gene TPM
      (TcgaTargetGtex_rsem_gene_tpm.gz, rows=Ensembl gene IDs, columns=samples)
      -- This file was actually downloaded by segmented_resume_download, not simulated data.

Workflow:
  1. Parse ZP3 and various immune marker gene symbols -> Ensembl IDs (obtained via Ensembl REST and cached to
     ensg_map.json, reproducible offline; fallback to built-in mapping table on failure).
  2. Stream-read the gz file and keep only target gene rows (avoid loading the full 1.3GB into memory).
  3. Filter TCGA tumor samples (sample ID pattern TCGA-XXXX-####-01) and group by cancer type.
  4. For each cancer type and each immune gene set, compute the
     Spearman rank correlation between ZP3 expression and the gene-set mean score (expression is log2(TPM) right-skewed, so use Spearman rather than Pearson).
  5. Generate a pan-cancer heatmap + association strength bar charts for each cancer type, and save detailed/summary CSVs.

Note: expression values are log2(TPM); "-9.9658" is the Xena underflow floor for 0, retained as low expression rank
    (Rank correlation is robust to monotonic transformations and does not involve absolute-value interpretation).
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
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(BASE)))  # Project root
DATA = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling",
                    "TcgaTargetGtex_rsem_gene_tpm.gz")
CACHE = os.path.join(BASE, "ensg_map.json")

# ---------------------------------------------------------------------------
# 1. Immune gene sets (symbol definitions, consistent with the previous version)
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

# Built-in fallback table (used only when Ensembl REST is unavailable; covers the most critical genes)
FALLBACK_ENSG = {
    "ZP3": "ENSG00000188372", "CD163": "ENSG00000177697", "CD68": "ENSG00000097273",
    "CD274": "ENSG00000120217", "PDCD1": "ENSG00000188389", "CTLA4": "ENSG00000163558",
    "FOXP3": "ENSG00000049768", "IFNG": "ENSG00000111537", "TGFB1": "ENSG00000105329",
    "IL10": "ENSG00000136634", "TIGIT": "ENSG00000180554", "LAG3": "ENSG00000179868",
    "VSIR": "ENSG00000116286", "BTLA": "ENSG00000120735", "IDO1": "ENSG00000131203",
    "TREM2": "ENSG00000183463",
}


def get_ensg_map(symbols):
    """symbol -> Ensembl gene id (version removed). Read cache first; otherwise query Ensembl REST and cache."""
    if os.path.exists(CACHE):
        with open(CACHE) as f:
            cache = json.load(f)
    else:
        cache = {}
    needed = [s for s in symbols if s not in cache]
    if needed:
        import urllib.request
        print(f"  Resolving {len(needed)} symbols -> ENSG from Ensembl REST ...")
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
                    print(f"    !! {s} failed to resolve: {e}")
            time.sleep(0.05)
        with open(CACHE, "w") as f:
            json.dump(cache, f, indent=2)
    return {s: cache.get(s) for s in symbols}


# ---------------------------------------------------------------------------
# 2. Stream-read target genes
# ---------------------------------------------------------------------------
def get_tcga_disease_map():
    """Map TCGA participant barcode -> cancer type (e.g. BRCA/LUAD/LGG/GBM).
    Sample IDs look like TCGA-<TSS>-<participant>-<sample>, where the second segment is the tissue source site (TSS),
    not the cancer type; the cancer type must be obtained from the GDC case->project mapping. The result is cached to
    tcga_disease_map.json, reproducible offline."""
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
            print(f"  !! GDC mapping fetch failed: {e}")
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
        time.sleep(0.05)  # avoid triggering GDC rate limiting
    with open(cache, "w") as f:
        json.dump(cases, f, indent=2)
    print(f"  GDC mapping cached {len(cases)} participants->cancer types")
    return cases


def read_target_genes(path, target_ids, chunk=2000):
    """Stream-read gz TPM matrix, keep only target_ids (Ensembl, with or without version).
    Returns DataFrame: index=Ensembl id, columns=sample id."""
    target_strip = {t.split(".")[0]: t for t in target_ids}
    rows = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        kept_samples = samples  # keep all for now, filter TCGA later
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
                print(f"    Scanned {n} gene lines...")
    df = pd.DataFrame.from_dict(rows, orient="index", columns=kept_samples)
    return df


# ---------------------------------------------------------------------------
# 3. Main pipeline
# ---------------------------------------------------------------------------
def main():
    print("=== TCGA pan-cancer analysis: tissue specificity of ZP3 and immune features (real data) ===\n")
    if not os.path.exists(DATA):
        print(f"!! Real data file not found: {DATA}\n    Please run segmented_resume_download to download "
              f"TcgaTargetGtex_rsem_gene_tpm.gz")
        sys.exit(1)

    all_symbols = [ZP3_SYMBOL] + sorted({g for s in IMMUNE_GENE_SETS.values() for g in s})
    sym2ensg = get_ensg_map(all_symbols)
    unresolved = [s for s, e in sym2ensg.items() if not e]
    if unresolved:
        print(f"  !! Unresolved genes (will be skipped): {unresolved}")
    # Target Ensembl id list
    target_ids = [e for e in sym2ensg.values() if e]
    print(f"  Parsed {len(target_ids)} Ensembl genes (including ZP3)\n")

    print("2. Streaming read real TPM matrix (target genes only)...")
    t0 = time.time()
    mat = read_target_genes(DATA, target_ids)
    print(f"  Reading complete, matrix {mat.shape[0]} genes × {mat.shape[1]} samples,"
          f"Elapsed time {time.time()-t0:.1f}s")

    # ZP3 row
    zp3_ensg = sym2ensg.get(ZP3_SYMBOL)
    if zp3_ensg not in mat.index:
        print(f"!! No ZP3 ({zp3_ensg}) in matrix, exiting"); sys.exit(1)
    zp3_vec_all = mat.loc[zp3_ensg]

    # Filter TCGA tumor samples
    samples = list(mat.columns)
    tcga_mask = [s.startswith("TCGA-") and s.split("-")[3].startswith("01") for s in samples]
    tcga_samples = [s for s, m in zip(samples, tcga_mask) if m]
    print(f"  TCGA tumor sample count: {len(tcga_samples)}")
    mat_t = mat[tcga_samples]

    # Group by cancer type: sample barcode = TCGA-<TSS>-<participant>-<sample>,
    # Field 2 is tissue source site (TSS), not cancer type; use GDC participant -> cancer type mapping.
    participant_of = {s: "-".join(s.split("-")[:3]) for s in tcga_samples}
    p2cancer = get_tcga_disease_map()  # participant barcode -> cancer type abbreviation
    cancer_of = {}
    for s in tcga_samples:
        cancer_of[s] = p2cancer.get(participant_of[s], "UNKNOWN")
    cancers = {}
    for s in tcga_samples:
        cancers.setdefault(cancer_of[s], []).append(s)
    cancers = {c: v for c, v in cancers.items() if c != "UNKNOWN"}
    n_unknown = sum(1 for c in cancer_of.values() if c == "UNKNOWN")
    if n_unknown:
        print(f"  (Warning: {n_unknown} samples could not be mapped to cancer types, excluded)")

    print(f"\n3. Calculate ZP3-immune association per cancer type (Spearman rho)...")
    print("   Scoring method: z-score consensus method (each gene standardized across samples then averaged)")
    records = []
    for cancer, sams in cancers.items():
        if len(sams) < 30:
            continue
        zp3_v = mat_t.loc[zp3_ensg, sams].values.astype(float)
        for set_name, syms in IMMUNE_GENE_SETS.items():
            ensgs = [sym2ensg[s] for s in syms if sym2ensg.get(s) in mat_t.index]
            if not ensgs:
                continue
            # z-score consensus method: standardize each gene across samples first, then average
            sub = mat_t.loc[ensgs, sams]           # genes x samples
            gene_mean = sub.mean(axis=1)           # mean per gene
            gene_std = sub.std(axis=1)             # per-gene standard deviation
            # Filter genes with std=0 (e.g., all are baseline values)
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
        print("!! No valid association results (check gene parsing/sample sizes)"); sys.exit(1)
    print(f"  Total {len(res)} (cancer type×gene set) associations, covering {res['Cancer_Code'].nunique()} cancer types")

    # Mean association per cancer type
    cancer_summary = res.groupby("Cancer_Code").agg(
        Avg_Rho=("Rho", "mean"),
        Sig_Count=("Significant", "sum"),
        N_sets=("Feature", "count"),
    ).reset_index().sort_values("Avg_Rho", ascending=False)
    cancer_summary["Sig_Count"] = cancer_summary["Sig_Count"].astype(int)

    print("\n  ZP3-immune association strength ranking (Top 8):")
    print("  " + "-" * 58)
    for _, row in cancer_summary.head(8).iterrows():
        print(f"  {row['Cancer_Code']:6s} | ρ = {row['Avg_Rho']:+.3f} | "
              f"Sig sets: {row['Sig_Count']}/{int(row['N_sets'])} | n={row['N_sets']*0}")

    # heatmap
    print("\n4. Generating pancancer association heatmap and bar chart...")
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
    print("  Saved fig_tcga_pancan_zp3_heatmap.png")

    res.to_csv(os.path.join(BASE, "tcga_pancan_zp3_correlations.csv"), index=False)
    cancer_summary.to_csv(os.path.join(BASE, "tcga_pancan_cancer_summary.csv"), index=False)
    print("  Saved tcga_pancan_zp3_correlations.csv / tcga_pancan_cancer_summary.csv")

    # conclusion
    print("\n5. Conclusion (real data):")
    strong = cancer_summary[cancer_summary["Avg_Rho"] > 0.15]["Cancer_Code"].tolist()
    moderate = cancer_summary[(cancer_summary["Avg_Rho"] > 0.05) &
                              (cancer_summary["Avg_Rho"] <= 0.15)]["Cancer_Code"].tolist()
    weak = cancer_summary[cancer_summary["Avg_Rho"] <= 0.05]["Cancer_Code"].tolist()
    print(f"  Strong correlation (ρ>0.15): {len(strong)} types — {', '.join(strong[:8])}")
    print(f"  Moderate correlation (0.05<ρ≤0.15): {len(moderate)} types — {', '.join(moderate[:8])}")
    print(f"  Weak/no correlation (ρ≤0.05): {len(weak)} types — {', '.join(weak[:8])}")
    print("\n=== Analysis complete (real TCGA TARGET GTEx data, replacing the original simulated data)===")


if __name__ == "__main__":
    main()
