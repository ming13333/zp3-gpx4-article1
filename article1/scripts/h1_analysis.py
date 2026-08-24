# -*- coding: utf-8 -*-
"""
H1 recon: survey of ZP3 expression in CNS immune cells (TAM/microglia/DC) in the glioma single-cell atlas
Data source: GSE141982 (TISCH GBM, ELSA classifier reference) — 10x mtx downloaded from GEO
Pipeline: load -> QC -> normalize/log1p -> HVG -> PCA -> neighbors -> UMAP -> Leiden
      -> score cell types against reference markers -> compute the fraction of ZP3+ cells
Output: table of ZP3 expression by cell type + UMAP plots (including ZP3 expression) + cell type annotation plot
"""
import os, warnings
import numpy as np
import pandas as pd
import scanpy as sc

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output", "h1_pilot")
warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
sc.settings.figdir = BASE
np.random.seed(0)

RAW = os.path.join(BASE, "raw")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article1", "results", "h1_pilot")

# ---------------------------------------------------------------
# 1. Load 4 samples (10x mtx)
# ---------------------------------------------------------------
# After tar extraction the files are flat; first file them into subdirectories by sample prefix
# Idempotent implementation: process only files (skip existing subdirectories), and already-filed subdirectories are
# skipped by the os.path.isfile guard — a second run won't crash when RAW is left with only subdirectories.
import re, shutil
prefixes = set()
for fn in os.listdir(RAW):
    fp = os.path.join(RAW, fn)
    if not os.path.isfile(fp):
        continue
    m = re.match(r"^(GSM\d+_[^_]+)_", fn)
    if m:
        prefixes.add(m.group(1))
for pre in prefixes:
    dstdir = os.path.join(RAW, pre)
    os.makedirs(dstdir, exist_ok=True)
    for fn in os.listdir(RAW):
        fp = os.path.join(RAW, fn)
        if not os.path.isfile(fp):
            continue
        if fn.startswith(pre + "_"):
            shutil.move(fp, os.path.join(dstdir, fn))

sample_dirs = sorted([d for d in os.listdir(RAW) if os.path.isdir(os.path.join(RAW, d))])
print("Found samples:", sample_dirs)
adatas = []
for d in sample_dirs:
    path = os.path.join(RAW, d)
    ad = sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=True, prefix=d+"_")
    ad.obs["sample"] = d
    adatas.append(ad)
    print(f"  {d}: {ad.n_obs} cells x {ad.n_vars} genes")

adata = sc.concat(adatas, join="outer", label="sample")
print("After merge:", adata.n_obs, "cells x", adata.n_vars, "genes")

# ---------------------------------------------------------------
# 2. QC
# ---------------------------------------------------------------
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
print("Cell count before QC:", adata.n_obs)
# Filter: at least 200 genes per cell; each gene expressed in at least 3 cells; mitochondrial < 25%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["pct_counts_mt"] < 25].copy()
print("Cell count after QC:", adata.n_obs)

# ---------------------------------------------------------------
# 3. Normalization / log / HVG / PCA / neighbors / UMAP / Leiden
# ---------------------------------------------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30, random_state=0)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0)
sc.tl.umap(adata, random_state=0)
sc.tl.leiden(adata, resolution=0.6, random_state=0, key_added="leiden")
print(" Leiden cluster count:", adata.obs["leiden"].nunique())

# ---------------------------------------------------------------
# 4. Cell type annotation via reference marker scoring
#    Use "cluster-level annotation" -- for each Leiden cluster, assign based on the
#    mean expression of lineage markers, avoiding noise contamination of per-cell max-score in microglia/myeloid.
# ---------------------------------------------------------------
markers = {
    "Microglia":      ["CX3CR1","P2RY12","TMEM119","SALL1","CSF1R","ITGAM","SPI1"],
    "TAM_Macrophage": ["CD68","CD163","LYZ","SPI1","C1QA","C1QB","CD14","VSIG4","MRC1"],
    "DC":             ["CLEC9A","ITGAX","HLA-DRA","FCER1A","LAMP3","CITED1","CD1C","BATF3"],
    "T_cell":         ["CD3D","CD3E","CD8A","CD4","IL7R","TRAC"],
    "B_cell":         ["CD79A","MS4A1","IGHM","CD19"],
    "NK":             ["NCAM1","NKG7","KLRD1","KLRF1"],
    "Malignant":      ["EGFR","PDGFRA","SOX2","OLIG2","CHI3L1","VIM","TOP2A","MKI67"],
    "Oligodendrocyte":["MBP","PLP1","OLIG1","MOBP"],
    "Astrocyte":      ["GFAP","AQP4","SLC1A2","GJA1","ALDH1L1"],
    "Endothelial":    ["PECAM1","CLDN5","VWF","FLT1","PLVAP"],
    "Pericyte":       ["PDGFRB","RGS5","ACTA2"],
}
valid_markers = {k: [g for g in v if g in adata.raw.var_names] for k, v in markers.items()}

# Per-cell scores (for UMAP display)
for k, v in valid_markers.items():
    if v:
        sc.tl.score_genes(adata, v, score_name=f"score_{k}")

# Cluster-level annotation: mean expression of each lineage marker per leiden cluster
raw_expr = adata.raw.X  # sparse (cells x genes, log-normalized)
gene_names = list(adata.raw.var_names)
marker_idx = {k: [gene_names.index(g) for g in v] for k, v in valid_markers.items() if v}
clusters = adata.obs["leiden"].unique()
cluster_map = {}
cluster_marker_mean = {}
for cl in clusters:
    cells = np.where(adata.obs["leiden"].values == cl)[0]
    sub = raw_expr[cells, :]
    means = {}
    for k, idx in marker_idx.items():
        means[k] = float(sub[:, idx].mean())
    best = max(means, key=means.get)
    cluster_map[cl] = best if means[best] > 0 else "Unassigned"
    cluster_marker_mean[cl] = means
adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_map).values

# Fix per-cell score column names (remove "/" for h5ad saving)
for k in valid_markers:
    if f"score_{k}" in adata.obs.columns:
        adata.obs.rename(columns={f"score_{k}": f"score_{k.replace('/','_')}"}, inplace=True)

# ---------------------------------------------------------------
# 5. ZP3 expression statistics (based on raw normalized data)
# ---------------------------------------------------------------
zp3 = "ZP3"
if zp3 in adata.raw.var_names:
    zp3_expr = adata.raw[:, zp3].X.toarray().ravel()  # log1p normalized counts
else:
    zp3_expr = np.zeros(adata.n_obs)
adata.obs["ZP3_expr"] = zp3_expr
adata.obs["ZP3_pos"] = zp3_expr > 0   # count>0 means detected

print("\n==================== H1 Results: ZP3 Expression Survey ====================")
print("Total cell count (post-QC):", adata.n_obs)
print("ZP3 overall positive rate: {:.2f}%".format(100*adata.obs['ZP3_pos'].mean()))

# ZP3 expression table by cluster-level annotation
rows = []
for ct in sorted(adata.obs["cell_type"].unique()):
    sub = adata.obs[adata.obs["cell_type"] == ct]
    n = len(sub)
    if n == 0: continue
    pos = sub["ZP3_pos"].sum()
    rows.append({
        "cell_type": ct, "n_cells": int(n),
        "pct_ZP3_positive": round(100*pos/n, 2),
        "mean_ZP3_logexpr": round(float(sub["ZP3_expr"].mean()), 4),
        "n_ZP3_positive": int(pos),
    })
res_df = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
print("\n--- ZP3 expression by annotated cell type (cluster level) ---")
print(res_df.to_string(index=False))

# Summary of central immune cells (TAM/Microglia/DC)
immune_types = ["TAM_Macrophage","Microglia","DC"]
imm = adata.obs[adata.obs["cell_type"].isin(immune_types)]
print("\n--- Central immune cell summary (TAM/Microglia/DC) ---")
print("Cell count:", len(imm), " ZP3+:", int(imm['ZP3_pos'].sum()),
      " Positive rate: {:.2f}%".format(100*imm['ZP3_pos'].mean()))

# ZP3 positive rate per sample (batch check)
print("\n--- ZP3 positive rate per sample ---")
for s in sorted(adata.obs["sample"].unique()):
    sub = adata.obs[adata.obs["sample"] == s]
    print(f"  {s}: n={len(sub)}, ZP3+={int(sub['ZP3_pos'].sum())}, {100*sub['ZP3_pos'].mean():.2f}%")

# Details within myeloid clusters
print("\n--- ZP3 in myeloid-related Leiden clusters (cluster-level check) ---")
for cl in sorted(clusters, key=lambda x: int(x)):
    if cluster_map[cl] in immune_types:
        sub = adata.obs[adata.obs["leiden"] == cl]
        print(f"  cluster {cl} -> {cluster_map[cl]}: n={len(sub)}, ZP3+={int(sub['ZP3_pos'].sum())}, {100*sub['ZP3_pos'].mean():.2f}%")

# ---------------------------------------------------------------
# 6. Visualization
# ---------------------------------------------------------------
sc.pl.umap(adata, color="cell_type", title="Cell types (GSE141982 GBM)",
           save="_celltypes.png", show=False, frameon=False, legend_fontsize=8)
sc.pl.umap(adata, color=["ZP3_expr","ZP3_pos"], title="ZP3 expression",
           save="_ZP3.png", show=False, frameon=False, cmap="viridis")
sc.pl.umap(adata, color="leiden", save="_leiden.png", show=False, frameon=False, legend_fontsize=8)

# Save
res_df.to_csv(os.path.join(OUT, "h1_zp3_by_celltype.csv"), index=False)
imm_summary = pd.DataFrame([{
    "group":"CNS_immune_TAM_MG_DC","n_cells":int(len(imm)),
    "n_ZP3_positive":int(imm['ZP3_pos'].sum()),
    "pct_ZP3_positive":round(100*imm['ZP3_pos'].mean(),2),
    "mean_ZP3_logexpr":round(float(imm['ZP3_expr'].mean()),4)}])
imm_summary.to_csv(os.path.join(OUT, "h1_zp3_immune_summary.csv"), index=False)
adata.write(os.path.join(OUT, "h1_adata.h5ad"))
print("\nSaved: h1_zp3_by_celltype.csv, h1_zp3_immune_summary.csv, h1_adata.h5ad, 3 UMAP plots")
