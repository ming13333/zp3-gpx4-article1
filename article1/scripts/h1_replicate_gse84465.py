# -*- coding: utf-8 -*-
"""
H1 独立复现（GSE84465 — GBM 单细胞，genes×cells 矩阵，空格分隔）
无注释列 -> 用泛髓系门控 + MG/TAM/DC marker 打分细分，重算 ZP3 阳性率。
目标：在第二个独立胶质瘤单细胞数据集确认 ZP3 在髓系（尤其 MG-like）富集。
"""
import os, gzip, numpy as np, pandas as pd

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article1", "results", "h1_pilot")
PATH = os.path.join(OUT, "GSE84465_GBM_All_data.csv.gz")

# 1. 探测分隔符（空格 vs 逗号）
with gzip.open(PATH, "rt") as f:
    first = f.readline()
sep = r"\s+" if first.count(" ") > first.count(",") else ","
print(f"分隔符: {sep!r}")

# 2. 读取（基因×细胞）
df = pd.read_csv(PATH, sep=sep, index_col=0, compression="gzip")
print(f"原始 shape (genes×cells): {df.shape}")
assert "ZP3" in df.index, "表达矩阵无 ZP3 行"

# 转置 -> 细胞×基因
expr = df.T
print(f"表达矩阵 (cells×genes): {expr.shape}")

# 3. ZP3 总体
zp3 = pd.to_numeric(expr["ZP3"], errors="coerce")
print(f"\n总体 ZP3 阳性率(>0): {100*(zp3>0).mean():.2f}% | 均值: {zp3.mean():.3f}")

# 4. marker 打分细分髓系
pan_myeloid = ["CD68", "LYZ", "C1QA", "C1QB", "ITGAM", "CSF1R", "CD14"]
MG  = ["CX3CR1", "P2RY12", "TMEM119", "SALL1", "SIGLEC11"]
TAM = ["CD163", "VSIG4", "MRC1", "MSR1", "FOLR2"]
DC  = ["CLEC9A", "FCER1A", "CD1C", "LAMP3", "BATF3", "ITGAX"]

def cmean(gl):
    ix = [g for g in gl if g in expr.columns]
    return expr[ix].mean(axis=1) if ix else pd.Series(0.0, index=expr.index)

pm = cmean(pan_myeloid); mg = cmean(MG); tam = cmean(TAM); dc = cmean(DC)
myeloid = pm > 0
print(f"\n泛髓系门控细胞数: {int(myeloid.sum())} ({100*myeloid.mean():.1f}%)")

M = pd.DataFrame({"MG": mg, "TAM": tam, "DC": dc})[myeloid]
best = M.idxmax(axis=1); val = M.max(axis=1)
subclass = best.where(val > 0, "Unassigned")

print("\n--- 细分髓系 ZP3 (GSE84465, 独立复现) ---")
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
print(f"\n髓系池总体: n={int(myeloid.sum())}, ZP3+={int((myeloid_zp3>0).sum())}, "
      f"{100*(myeloid_zp3>0).mean():.2f}%")
print(f"全数据集总体 ZP3+ 率: {100*(zp3>0).mean():.2f}%")
print("=== GSE84465 H1 复现完成 ===")
