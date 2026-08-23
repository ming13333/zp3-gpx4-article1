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
    return _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
ROOT = _project_root()
import os
# -*- coding: utf-8 -*-
"""
检验B: bulk GBM/LGG 中 ZP3↔TREM2 关联是否独立于"总髓系负荷"
- 用 GET 端点拉泛髓系 marker (CD68,CD14,LYZ,CSF1R,ITGAM) 到 expr
- 髓系指数 = 可用 marker 的 z-scored mean
- 粗相关 vs 偏相关(控制髓系指数)
GET 端点（样本列表用 *_rna_seq_v2_mrna），逐基因。
"""
import requests, os, time
import numpy as np, pandas as pd, scipy.stats as st
API="https://www.cbioportal.org/api"
BASE=os.path.join(ROOT, "output", "h2_bulk")
H={"Accept":"application/json"}

def get(url,params=None,retries=3):
    for i in range(retries):
        try:
            r=requests.get(url,params=params,headers=H,timeout=120);r.raise_for_status();return r.json()
        except Exception as e:
            if i==retries-1: print("ERR",e);return []
            time.sleep(2)

def entrez(sym):
    d=get("%s/genes/%s"%(API,sym))
    return d.get("entrezGeneId") if isinstance(d,dict) else None

marker_syms=["CD68","CD14","LYZ","CSF1R","ITGAM"]
for study,suffix in [("gbm_tcga","gbm"),("lgg_tcga","lgg")]:
    mp="%s_rna_seq_v2_mrna"%study
    sl="%s_rna_seq_v2_mrna"%study
    expr=pd.read_csv(os.path.join(BASE,"expr_%s_tcga_patient.csv"%suffix),index_col=0)
    got=[]
    for sym in marker_syms:
        e=entrez(sym)
        if not e: continue
        zs=get("%s/molecular-profiles/%s/molecular-data"%(API,mp),params={"sampleListId":sl,"entrezGeneId":e})
        n=0
        for r in zs:
            pid=r.get("patientId");val=r.get("value")
            if pid and val is not None:
                expr.loc[pid,sym]=float(val); n+=1
        if n>30: got.append(sym)
        time.sleep(0.2)
    print("=== %s: 拉到的髓系 marker %s ==="%(study,got))
    if got:
        zs_=(expr[got]-expr[got].mean())/expr[got].std()
        expr["myeloid_idx"]=zs_.mean(axis=1)
        sub=expr.dropna(subset=["ZP3","TREM2","myeloid_idx"])
        r0,p0=st.pearsonr(sub["ZP3"],sub["TREM2"])
        r1,p1=st.pearsonr(sub["ZP3"],sub["myeloid_idx"])
        r2,p2=st.pearsonr(sub["TREM2"],sub["myeloid_idx"])
        def reg(y,x):
            X=np.column_stack([np.ones(len(x)),x]);c,_,_,_=np.linalg.lstsq(X,y,rcond=None);return c[0]+X[:,1]*c[1]
        rzp=sub["ZP3"].values-reg(sub["ZP3"].values,sub["myeloid_idx"].values)
        rtr=sub["TREM2"].values-reg(sub["TREM2"].values,sub["myeloid_idx"].values)
        rp,pp=st.pearsonr(rzp,rtr)
        print("  n=%d | 髓系指数基因: %s"%(len(sub),got))
        print("  ZP3 vs 髓系指数:      r=%.3f p=%.3g"%(r1,p1))
        print("  TREM2 vs 髓系指数:    r=%.3f p=%.3g"%(r2,p2))
        print("  ZP3 vs TREM2 粗相关:  r=%.3f p=%.3g"%(r0,p0))
        print("  ZP3 vs TREM2 偏相关(控制髓系指数): r_partial=%.3f p=%.3g"%(rp,pp))
        print("  => 控制髓系负荷后 %s"%(("仍显著: 独立于总髓系量的特异关联" if pp<0.05 else "不显著: 表观相关主要由总髓系负荷驱动(混杂)")))
        expr.to_csv(os.path.join(BASE,"expr_%s_b_myeloid.csv"%suffix))
