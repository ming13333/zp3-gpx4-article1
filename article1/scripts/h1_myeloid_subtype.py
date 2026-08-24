# -*- coding: utf-8 -*-
"""
H1 Deep-dive: subdivide myeloid and recalculate ZP3
- Load h1_adata.h5ad (GSE141982 GBM, already QC)
- Pan-myeloid gating defines the myeloid pool, avoiding missing myeloid cells hidden in NK clusters if only TAM clusters are used
- Use scores for homeostatic microglia (MG) vs BMDM/TAM vs DC markers to subdivide subtypes
- Recalculate ZP3+ rate separately to determine enrichment location
"""
import os, numpy as np, pandas as pd, scanpy as sc

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article1", "results", "h1_pilot")
ad = sc.read_h5ad(os.path.join(OUT, "h1_adata.h5ad"))
ad.obs_names_make_unique()

raw = ad.raw[:, ad.raw.var_names].X  # sparse, log1p(normalized)
genes = list(ad.raw.var_names)
def gidx(glist):
    return [genes.index(g) for g in glist if g in genes]

pan_myeloid = ["CD68", "LYZ", "C1QA", "C1QB", "ITGAM", "CSF1R", "CD14"]
MG  = ["CX3CR1", "P2RY12", "TMEM119", "SALL1", "SIGLEC11"]
TAM = ["CD163", "VSIG4", "MRC1", "MSR1", "FOLR2"]
DC  = ["CLEC9A", "FCER1A", "CD1C", "LAMP3", "BATF3", "ITGAX"]

def cellmean(glist):
    ix = gidx(glist)
    return np.asarray(raw[:, ix].mean(axis=1)).ravel()

ad.obs["pan_myeloid"] = cellmean(pan_myeloid)
ad.obs["MG_score"]    = cellmean(MG)
ad.obs["TAM_score"]   = cellmean(TAM)
ad.obs["DC_score"]    = cellmean(DC)

# ZP3
zp3 = ad.raw[:, "ZP3"].X.toarray().ravel()
ad.obs["ZP3_expr"] = zp3
ad.obs["ZP3_pos"]  = zp3 > 0

# ---- 1) Diagnostics: mean scores of myeloid/subtype per leiden cluster, to see if myeloid cells are hidden in NK clusters ----
print("="*72)
print("Diagnostics: mean myeloid marker per leiden cluster (log1p)")
diag = ad.obs.groupby("leiden").agg(
    n=("ZP3_pos", "size"),
    pan_myeloid=("pan_myeloid", "mean"),
    MG=("MG_score", "mean"),
    TAM=("TAM_score", "mean"),
    DC=("DC_score", "mean"),
    ZP3_pos_rate=("ZP3_pos", "mean"),
).round(3)
print(diag.to_string())

# ---- 2) Pan-myeloid gating defines the myeloid pool ----
# Gating: pan_myeloid > 0 (at least one pan-myeloid marker detectable)
myeloid = ad.obs["pan_myeloid"] > 0
ad.obs["myeloid"] = myeloid
print("\nMyeloid pool cell count after pan-myeloid gating:", int(myeloid.sum()),
      f"({100*myeloid.mean():.1f}% of {ad.n_obs})")

# ---- 3) Subclassification within the myeloid pool ----
sub = ad.obs[myeloid]
M = sub[["MG_score", "TAM_score", "DC_score"]].values
best = np.argmax(M, axis=1)
val = M[np.arange(len(M)), best]
classes = np.array(["MG", "TAM", "DC"])[best]
classes[val <= 0] = "Unassigned"
ad.obs.loc[myeloid, "myeloid_subclass"] = classes

# ---- 4) Recalculate ZP3+ rates separately ----
print("\n" + "="*72)
print("Myeloid subclass ZP3 expression (GSE141982 GBM)")
rows = []
for cls in ["MG", "TAM", "DC", "Unassigned"]:
    cc = ad.obs[myeloid & (ad.obs["myeloid_subclass"] == cls)]
    n = len(cc)
    if n == 0:
        continue
    pos = int(cc["ZP3_pos"].sum())
    rows.append({
        "myeloid_subclass": cls,
        "n_cells": n,
        "n_ZP3_pos": pos,
        "pct_ZP3_pos": round(100*pos/n, 2),
        "mean_ZP3_logexpr": round(float(cc["ZP3_expr"].mean()), 4),
    })
res = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
print(res.to_string(index=False))
print(f"\nMyeloid pool overall: n={int(myeloid.sum())}, ZP3+={int(ad.obs[myeloid]['ZP3_pos'].sum())}, "
      f"{100*ad.obs[myeloid]['ZP3_pos'].mean():.2f}%")
print(f"Full dataset overall ZP3+ rate: {100*ad.obs['ZP3_pos'].mean():.2f}%")

# ---- 5) Comparison with original cluster-level TAM annotation ----
orig_tam = ad.obs[ad.obs["cell_type"] == "TAM_Macrophage"]
print(f"\nOriginal cluster-level TAM cluster (leiden9): n={len(orig_tam)}, ZP3+={orig_tam['ZP3_pos'].sum()}, "
      f"{100*orig_tam['ZP3_pos'].mean():.2f}%")

# ---- 6) Save ----
res.to_csv(os.path.join(OUT, "h1_zp3_myeloid_subtype.csv"), index=False)
# UMAP coloring
sc.pl.umap(ad, color=["myeloid_subclass", "ZP3_expr"],
           save="_myeloid_subtype.png", show=False, frameon=False, legend_fontsize=8, cmap="viridis")
ad.write(os.path.join(OUT, "h1_adata_subtyped.h5ad"))
print("\nSaved: h1_zp3_myeloid_subtype.csv, h1_adata_subtyped.h5ad, umap_myeloid_subtype.png")
