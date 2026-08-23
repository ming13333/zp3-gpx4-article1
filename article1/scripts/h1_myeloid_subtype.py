# -*- coding: utf-8 -*-
"""
H1 深化：细分髓系重算 ZP3
- 读 h1_adata.h5ad (GSE141982 GBM, 已 QC)
- 泛髓系门控 (pan-myeloid) 定义髓系池，避免仅看 TAM 簇漏掉混入 NK 簇的髓系
- 用稳态小胶质 (MG) vs BMDM/TAM vs DC 标记打分，细分亚类
- 分别重算 ZP3+ 率，判断富集落点
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

# ---- 1) 诊断：每个 leiden 簇的髓系/亚类平均得分，看髓系是否藏在 NK 簇 ----
print("="*72)
print("诊断：各 leiden 簇的髓系标记均值 (log1p)")
diag = ad.obs.groupby("leiden").agg(
    n=("ZP3_pos", "size"),
    pan_myeloid=("pan_myeloid", "mean"),
    MG=("MG_score", "mean"),
    TAM=("TAM_score", "mean"),
    DC=("DC_score", "mean"),
    ZP3_pos_rate=("ZP3_pos", "mean"),
).round(3)
print(diag.to_string())

# ---- 2) 泛髓系门控定义髓系池 ----
# 门控：pan_myeloid > 0（至少一个泛髓系标记可检测）
myeloid = ad.obs["pan_myeloid"] > 0
ad.obs["myeloid"] = myeloid
print("\n泛髓系门控后髓系池细胞数:", int(myeloid.sum()),
      f"({100*myeloid.mean():.1f}% of {ad.n_obs})")

# ---- 3) 髓系池内细分亚类 ----
sub = ad.obs[myeloid]
M = sub[["MG_score", "TAM_score", "DC_score"]].values
best = np.argmax(M, axis=1)
val = M[np.arange(len(M)), best]
classes = np.array(["MG", "TAM", "DC"])[best]
classes[val <= 0] = "Unassigned"
ad.obs.loc[myeloid, "myeloid_subclass"] = classes

# ---- 4) 分别重算 ZP3+ 率 ----
print("\n" + "="*72)
print("细分髓系 ZP3 表达（GSE141982 GBM）")
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
print(f"\n髓系池总体: n={int(myeloid.sum())}, ZP3+={int(ad.obs[myeloid]['ZP3_pos'].sum())}, "
      f"{100*ad.obs[myeloid]['ZP3_pos'].mean():.2f}%")
print(f"全数据集总体 ZP3+ 率: {100*ad.obs['ZP3_pos'].mean():.2f}%")

# ---- 5) 与原始簇水平 TAM 注释对比 ----
orig_tam = ad.obs[ad.obs["cell_type"] == "TAM_Macrophage"]
print(f"\n原始簇水平 TAM 簇(leiden9): n={len(orig_tam)}, ZP3+={orig_tam['ZP3_pos'].sum()}, "
      f"{100*orig_tam['ZP3_pos'].mean():.2f}%")

# ---- 6) 保存 ----
res.to_csv(os.path.join(OUT, "h1_zp3_myeloid_subtype.csv"), index=False)
# UMAP 着色
sc.pl.umap(ad, color=["myeloid_subclass", "ZP3_expr"],
           save="_myeloid_subtype.png", show=False, frameon=False, legend_fontsize=8, cmap="viridis")
ad.write(os.path.join(OUT, "h1_adata_subtyped.h5ad"))
print("\n已保存: h1_zp3_myeloid_subtype.csv, h1_adata_subtyped.h5ad, umap_myeloid_subtype.png")
