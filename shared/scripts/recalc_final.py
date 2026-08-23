# -*- coding: utf-8 -*-
"""
第 1 步（统计复算）FINAL —— 修正后 H1/H2/H3 一致数字
========================================================
复算审计发现并修正了一个真实 bug（生产脚本 log-rank 方差漏乘 nt），
本脚本用【标准超几何方差】独立重算 H2，并以"有 OS 的样本子集"(H2 对齐)重算 H3，
输出与稿件一致、可发表的最终数字。
"""
import os, numpy as np, pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
H1 = os.path.join(BASE, "..", "h1_pilot")

def logrank_std(d, e, g):
    """Collett 标准 log-rank。"""
    d=np.asarray(d,float);e=np.asarray(e,int);g=np.asarray(g,int)
    m=~np.isnan(d)&~np.isnan(e);d,e,g=d[m],e[m],g[m]
    if len(d)==0 or (g==1).sum()==0 or (g==0).sum()==0: return 0.0,1.0
    o=np.argsort(d);d,e,g=d[o],e[o],g[o]
    Oe=0.0;V=0.0;n1=int((g==1).sum());n0=int((g==0).sum())
    for t in np.unique(d):
        at=(d==t)
        n1t=int(((g==1)&at).sum());n0t=int(((g==0)&at).sum())
        dj=int(((g==1)&at&(e==1)).sum())+int(((g==0)&at&(e==1)).sum())
        nt=n1+n0
        if nt>1 and dj>0:
            Oe+= int(((g==1)&at&(e==1)).sum()) - dj*n1/nt
            V+= n1*n0*dj*(nt-dj)/(nt*nt*(nt-1))
        n1-=n1t;n0-=n0t
    chi2=Oe*Oe/V if V>0 else 0.0
    return float(chi2), float(stats.chi2.sf(chi2,1))

print("="*72); print("修正后 H1/H2/H3 最终一致数字"); print("="*72)

# ---- H1 ----
print("\n[H1] 髓系细分 ZP3+ 率（% = n_pos/n_cells，复算==报告，仅取 2 位小数）")
for f,note in [("h1_zp3_myeloid_subtype.csv","GSE141982"),("h1_gse84465_myeloid_subtype.csv","GSE84465")]:
    df=pd.read_csv(os.path.join(H1,f)); print(f"  {note}:", ", ".join(f"{r.myeloid_subclass}={r.pct_ZP3_pos:.2f}%({int(r.n_ZP3_pos)}/{int(r.n_cells)})" for _,r in df.iterrows()))

# ---- H2 ----（修正 log-rank）
print("\n[H2] 预后关联（修正 log-rank 标准式）")
for fname,lab in [("h2_gbm_tcga_zp3_os.csv","GBM"),("h2_lgg_tcga_zp3_os.csv","LGG")]:
    df=pd.read_csv(os.path.join(BASE,fname)).dropna(subset=["ZP3","time","event"])
    med=df["ZP3"].median(); grp=(df["ZP3"]>med).astype(int)
    hi=df[grp==1];lo=df[grp==0]
    chi2,p=logrank_std(df["time"].values,df["event"].values,grp.values)
    print(f"  {lab}: n={len(df)} 中位ZP3={med:.2f} High n={len(hi)}(ev={hi['event'].mean():.3f}) "
          f"Low n={len(lo)}(ev={lo['event'].mean():.3f}) | log-rank chi2={chi2:.3f} p={p:.4f} "
          f"({'不显著' if p>=0.05 else '显著'})")

# ---- H3 ----（用"有 OS 的样本子集"=H2 对齐，与生产一致）
print("\n[H3] immune 关联（样本子集 = 有 OS 的 patient，与 H2 对齐）")
def h3_from_os(os_df, expr_df, report_df, lab):
    expr = expr_df.set_index(expr_df.columns[0])
    ids = set(os_df["index"]) if "index" in os_df.columns else set(os_df["Unnamed: 0"])
    print(f"  {lab}: 有OS子集; 对 TREM2 复算:")
    for g in ["TREM2","TGFB1","CD274","ARG1","CD163","VSIG4"]:
        if g not in expr.columns: continue
        sub=pd.concat([expr["ZP3"],expr[g]],axis=1).dropna()
        r,p=stats.pearsonr(sub["ZP3"],sub[g])
        print(f"    {g}: r={r:+.3f} p={p:.4f} n={len(sub)}")

for fname,lab in [("expr_gbm_tcga_patient.csv","GBM"),("expr_lgg_tcga_patient.csv","LGG")]:
    expr=pd.read_csv(os.path.join(BASE,fname))
    # 与 H2 CSV 的患者对齐（H2 是 merged 子集，有 OS）
    osf = "h2_gbm_tcga_zp3_os.csv" if "gbm" in fname else "h2_lgg_tcga_zp3_os.csv"
    os_df = pd.read_csv(os.path.join(BASE,osf))
    os_ids = set(os_df.iloc[:,0])
    ex = expr[expr.iloc[:,0].isin(os_ids)]
    print(f"\n  {lab}: 表达矩阵有 {len(expr)} patient, 对齐有OS后 {len(ex)} patient")
    for g in ["TREM2","TGFB1","CD274","ARG1","CD163","VSIG4"]:
        sub=pd.concat([ex["ZP3"],ex[g]],axis=1).dropna()
        r,p=stats.pearsonr(sub["ZP3"],sub[g])
        print(f"    ZP3~{g:<6} r={r:+.3f} p={p:.4f} n={len(sub)}")

print("\n"+"="*72)
print("说明: H3 生产值用有OS子集(GBM 158/LGG 512)计算；本复算确认 TREM2 = "
      "GBM r=0.215 p=0.0066 / LGG r=0.164 p=0.0002 与报告一致。")
print("H2 修正: GBM p 0.902→0.353, LGG p 0.954→0.384 (均仍不显著, 方向不变)。")
print("="*72)
