#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
#37 Improvement ②: GPX4–ZP3 expression association (tying back to Cell 2026 GPX4-ZP3 immunosuppressive axis migration logic)
=====================================================================
Background: All three articles claim to migrate the Cell 2026 'extracellular GPX4–ZP3 immunosuppressive DAMP axis' from
reproductive biology to neuro-oncology, but the expression association between GPX4 and ZP3 has never been validated.
This script uses real data to fill in this back-reference evidence.

Data: local real TPM (TcgaTargetGtex_rsem_gene_tpm.gz, 1.3GB, already downloaded)
Target genes:
  - ZP3 (ENSG00000188372)
  - GPX4 (glutathione peroxidase 4, ENSG00000112715)  ← axis core
  - Ferroptosis scoring genes: GPX4 / ACSL4 / SLC7A11 / TFRC / FTL / FTH1 / NFE2L2 (NRF2)
    The score is constructed from two aspects: anti-ferroptosis (defense) and pro-ferroptosis (sensitivity)

Analysis content:
  1. GPX4–ZP3 expression Spearman correlation (all TCGA tumor samples + glioma GBM/LGG)
  2. Ferroptosis defense score (GPX4/SLC7A11/NFE2L2) association with ZP3
  3. Ferroptosis sensitivity score (ACSL4/TFRC) association with ZP3
  4. Pan-cancer stratification of GPX4–ZP3 association by cancer type, validating tissue specificity
  5. Output association plots + CSV, linking back to migration logic

Outputs: gpx4_zp3_*csv + fig_gpx4_zp3_association.png
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
# Target genes (symbol -> ensg)
# ---------------------------------------------------------------------------
GENES = {
    "ZP3": "ENSG00000188372",
    # Ferroptosis axis
    "GPX4": "ENSG00000112715",
    "ACSL4": "ENSG00000068366",
    "SLC7A11": "ENSG00000151012",
    "TFRC": "ENSG00000072274",      # iron uptake promotes ferroptosis
    "FTL": "ENSG00000087086",       # ferritin light chain (stores iron, protects against ferroptosis)
    "FTH1": "ENSG00000167996",      # ferritin heavy chain
    "NFE2L2": "ENSG00000116044",    # NRF2, ferroptosis defense transcription factor
    # downstream lipid peroxidation
    "ALOX15": "ENSG00000161905",
}

# Ferroptosis score decomposition
ANTI_FERROPTOSIS = ["GPX4", "SLC7A11", "NFE2L2", "FTL"]   # defense
PRO_FERROPTOSIS = ["ACSL4", "TFRC"]                       # sensitive/promotes death


def ensure_ensg(dict_path, extra):
    """Merge missing gene mappings into the cache and write them out."""
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
        print(f"  ensg_map.json supplemented with {sorted(added)}")
    return m


def read_target_genes(path, want_ensg, chunk=2000):
    """Stream-read gz TPM, keeping only target gene rows. Returns DataFrame: index=ensg, cols=sample."""
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
    print("=== #37 GPX4–ZP3 expression association analysis (real TPM data) ===\n")
    if not os.path.exists(DATA):
        print(f"!! missing TPM: {DATA}"); return

    # 1. Supplement ensg_map
    print("1. Ensure gene mapping (supplement ferroptosis genes like GPX4/ACSL4)...")
    ensg_map = ensure_ensg(ENSG_MAP, GENES)
    want_ensg = list(GENES.values())
    inv = {v: k for k, v in GENES.items()}

    # 2. Read target genes
    print("2. Stream-read actual TPM (target genes)...")
    t0 = time.time()
    mat = read_target_genes(DATA, want_ensg)
    print(f"    Matrix {mat.shape[0]} genes x {mat.shape[1]} samples, elapsed {time.time()-t0:.1f}s")
    missing = [sym for sym, e in GENES.items() if e not in mat.index]
    if missing:
        print(f"   !! Missing gene rows: {missing}")
    mat.index = [inv.get(e, e) for e in mat.index]

    # 3. TCGA tumor samples
    samples = list(mat.columns)
    tcga = [s for s in samples if s.startswith("TCGA-") and s.split("-")[3].startswith("01")]
    print(f"   TCGA tumor samples: {len(tcga)}")
    mat_t = mat[tcga]

    # 4. Cancer type mapping (mat_t is genes x samples; transpose and add cancer column)
    with open(DISEASE_MAP) as f:
        disease = json.load(f)
    def cancer_of(s):
        return disease.get("-".join(s.split("-")[:3]), "UNKNOWN")
    mat_tt = mat_t.T.copy()                # samples x genes
    mat_tt["cancer"] = [cancer_of(s) for s in mat_tt.index]
    mat_tt = mat_tt[mat_tt["cancer"] != "UNKNOWN"]

    zp3 = mat_tt["ZP3"]
    gpx4 = mat_tt["GPX4"]

    # ---- 5. Global + glioma GPX4–ZP3 ----
    print("\n3. GPX4–ZP3 association:")
    rec_global = []
    for label, sel in [
        ("ALL_TCGA", np.ones(len(mat_tt), dtype=bool)),
        ("GBM", mat_tt["cancer"].values == "GBM"),
        ("LGG", mat_tt["cancer"].values == "LGG"),
    ]:
        rho, p, n = spearman_p(zp3.values[sel], gpx4.values[sel])
        rec_global.append({"Cohort": label, "Rho": rho, "P": p, "N": n})
        print(f"   {label:10s}: GPX4–ZP3 ρ={rho:+.3f}, p={p:.3g}, n={n}")

    # ---- 6. Ferroptosis score with ZP3 ----
    def score_zs(genes, df):
        """z-score consensus: standardize each gene (column) across samples, then average genes per sample."""
        sub = df[genes].astype(float).T        # gene rows × sample columns
        valid = sub.std(axis=1) > 0            # drop std=0 genes
        if not valid.any():
            return pd.Series(np.nan, index=df.index)
        z = ((sub.loc[valid] - sub.loc[valid].mean(axis=1).values[:, None])
             / sub.loc[valid].std(axis=1).values[:, None])
        return z.mean(axis=0)                  # one score per sample

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

    # ---- 7. Pan-cancer GPX4–ZP3 stratification by cancer type ----
    rec_pan = []
    for c, grp in mat_tt.groupby("cancer"):
        if len(grp) < 30:
            continue
        rho, p, n = spearman_p(grp["ZP3"].values, grp["GPX4"].values)
        rec_pan.append({"Cancer": c, "Rho": rho, "P": p, "N": n})
    pan = pd.DataFrame(rec_pan).sort_values("Rho", ascending=False)
    print("\n4. Pan-cancer GPX4–ZP3 association (Top 10):")
    for _, r in pan.head(10).iterrows():
        star = "*" if r["P"] < 0.05 else " "
        print(f"   {r['Cancer']:6s}  ρ={r['Rho']:+.3f}  p={r['P']:.3g}  n={r['N']} {star}")

    # ---- 8. Save CSV ----
    pd.DataFrame(rec_global).to_csv(os.path.join(BASE, "gpx4_zp3_global.csv"), index=False)
    pd.DataFrame(rec_fx).to_csv(os.path.join(BASE, "gpx4_zp3_ferroptosis_score.csv"), index=False)
    pan.to_csv(os.path.join(BASE, "gpx4_zp3_pancancer.csv"), index=False)
    mat_tt[["ZP3", "GPX4", "cancer"]].to_csv(
        os.path.join(BASE, "gpx4_zp3_expr_matrix.csv"))

    # ---- 9. Plot ----
    print("\n5. Generating association plot...")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))

    # a) All TCGA scatter
    ax = axes[0]
    ax.scatter(gpx4.values, zp3.values, s=8, alpha=0.4, c="#378ADD")
    rho_all, p_all, _ = spearman_p(zp3.values, gpx4.values)
    ax.set_xlabel("GPX4 log2(TPM)")
    ax.set_ylabel("ZP3 log2(TPM)")
    ax.set_title(f"All TCGA (n={len(mat_tt)})\nρ={rho_all:+.3f}, p={p_all:.1e}", fontsize=11)

    # b) GBM/LGG scatter
    ax = axes[1]
    for c, col in [("GBM", "#A32D2D"), ("LGG", "#1D9E75")]:
        g = mat_tt[mat_tt["cancer"] == c]
        ax.scatter(g["GPX4"].values, g["ZP3"].values, s=12, alpha=0.6,
                   label=c, c=col)
    ax.set_xlabel("GPX4 log2(TPM)"); ax.set_ylabel("ZP3 log2(TPM)")
    ax.set_title("Glioma: GPX4–ZP3", fontsize=11)
    ax.legend()

    # c) Pan-cancer bar
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
    print("    Saved fig_gpx4_zp3_association.png")

    # ---- 10. Rebate migration logic conclusion ----
    print("\n=== Conclusion (referring back to Cell 2026 GPX4-ZP3 axis migration logic) ===")
    gbm_r = next(r for r in rec_global if r["Cohort"] == "GBM")
    lgg_r = next(r for r in rec_global if r["Cohort"] == "LGG")
    all_r = next(r for r in rec_global if r["Cohort"] == "ALL_TCGA")
    print(f"  · Global TCGA: GPX4–ZP3 ρ={all_r['Rho']:+.3f} (p={all_r['P']:.2g})")
    print(f"  · GBM: ρ={gbm_r['Rho']:+.3f} (p={gbm_r['P']:.2g})")
    print(f"  · LGG: ρ={lgg_r['Rho']:+.3f} (p={lgg_r['P']:.2g})")
    for r in rec_fx:
        print(f"  · ZP3–{r['Score']}: ρ={r['Rho']:+.3f} (p={r['P']:.2g})")
    print("\nConclusion recap: GPX4, as the core of the ferroptosis axis, if significantly positively correlated with ZP3,"
          "then supports『ferroptosis(ferroptosis) releases GPX4→secreted GPX4-ZP3 complex→"
          "the extension of the 'myeloid immunosuppression' Cell 2026 axis in neuro-oncology."
          "If GBM/LGG is consistent with the pan-cancer direction, it provides expression-level evidence for the migration logic.")
    print("\n=== #37 done ===")


if __name__ == "__main__":
    main()
