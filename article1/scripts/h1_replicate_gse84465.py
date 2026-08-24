# -*- coding: utf-8 -*-
"""
H1 independent replication (GSE84465 — GBM single-cell, genes×cells matrix, space-separated)
No annotation column -> use pan-myeloid gating + MG/TAM/DC marker scoring to subdivide, recalculate ZP3 positivity rate.
Goal: confirm ZP3 enrichment in myeloid cells (especially MG-like) in a second independent glioma single-cell dataset.
"""
import os, gzip, numpy as np, pandas as pd

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article1", "results", "h1_pilot")
PATH = os.path.join(OUT, "GSE84465_GBM_All_data.csv.gz")

# 1. Detect delimiter (space vs comma)
with gzip.open(PATH, "rt") as f:
    first = f.readline()
sep = r"\s+" if first.count(" ") > first.count(",") else ","
print(f"separator: {sep!r}")

# 2. Read (genes×cells)
df = pd.read_csv(PATH, sep=sep, index_col=0, compression="gzip")
print(f"Original shape (genes×cells): {df.shape}")
assert "ZP3" in df.index, "Expression matrix lacks ZP3 row"

# Transpose -> cells×genes
expr = df.T
print(f"Expression matrix (cells×genes): {expr.shape}")

# 3. ZP3 overall
zp3 = pd.to_numeric(expr["ZP3"], errors="coerce")
print(f"\nOverall ZP3 positivity rate (>0): {100*(zp3>0).mean():.2f}% | mean: {zp3.mean():.3f}")

# 4. marker scoring to subdivide myeloid
pan_myeloid = ["CD68", "LYZ", "C1QA", "C1QB", "ITGAM", "CSF1R", "CD14"]
MG  = ["CX3CR1", "P2RY12", "TMEM119", "SALL1", "SIGLEC11"]
TAM = ["CD163", "VSIG4", "MRC1", "MSR1", "FOLR2"]
DC  = ["CLEC9A", "FCER1A", "CD1C", "LAMP3", "BATF3", "ITGAX"]

def cmean(gl):
    ix = [g for g in gl if g in expr.columns]
    return expr[ix].mean(axis=1) if ix else pd.Series(0.0, index=expr.index)

pm = cmean(pan_myeloid); mg = cmean(MG); tam = cmean(TAM); dc = cmean(DC)
myeloid = pm > 0
print(f"\nPan-myeloid gated cell count: {int(myeloid.sum())} ({100*myeloid.mean():.1f}%)")

M = pd.DataFrame({"MG": mg, "TAM": tam, "DC": dc})[myeloid]
best = M.idxmax(axis=1); val = M.max(axis=1)
subclass = best.where(val > 0, "Unassigned")

print("\n--- Subcluster myeloid ZP3 (GSE84465, independent replication) ---")
rows = []
for cls in ["MG", "TAM", "DC", "Unassigned"]:
    cc = zp3[myeloid & (subclass == cls)]
    n = len(cc)
    if n == 0:
        continue
    rows.append({
        "myeloid_subclass": cls, "n_cells": n,
        "n_ZP3_pos": int((cc > 0).sum()),
        "pct_ZP3_pos": round(100*(cc > 0).mean(), 2),
        "mean_ZP3": round(float(cc.mean()), 4),
    })
res = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
print(res.to_string(index=False))
res.to_csv(os.path.join(OUT, "h1_gse84465_myeloid_subtype.csv"), index=False)

myeloid_zp3 = zp3[myeloid]
print(f"\nMyeloid pool overall: n={int(myeloid.sum())}, ZP3+={int((myeloid_zp3>0).sum())}, "
      f"{100*(myeloid_zp3>0).mean():.2f}%")
print(f"Full dataset overall ZP3+ rate: {100*(zp3>0).mean():.2f}%")
print("=== GSE84465 H1 reproduction complete ===")
