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
Method 3 + Test B fetch: CCLE pure tumor cell line ZP3 expression + bulk pan-myeloid index partial correlation
"""
import requests, time, sys, os, json
import numpy as np, pandas as pd

API = "https://www.cbioportal.org/api"
HEADERS = {"Accept":"application/json","Content-Type":"application/json"}
BASE = os.path.join(ROOT, "output", "h2_bulk")

def get(path, params=None, retries=3):
    for i in range(retries):
        try:
            r = requests.get(API+path, params=params, headers=HEADERS, timeout=45)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if i==retries-1:
                return {"error":str(e)}
            time.sleep(1.5)

def entrez(sym):
    d = get("/genes/%s"%sym)
    if isinstance(d,dict) and "entrezGeneId" in d: return d["entrezGeneId"]
    return None

# ---------------- Method 3: CCLE ----------------
print("="*70)
print("[Method3] CCLE pure tumor cell line ZP3 expression (test: whether pure tumor cells express ZP3)")
mp = "ccle_broad_2019_rna_seq_mrna"
# 1) sample list (all CCLE cell lines)
sl = get("/sample-lists/ccle_broad_2019_all/sample-ids")
if isinstance(sl, list) and sl:
    samples = sl
    print("  CCLE all cell lines:", len(samples))
    # 2) Pull ZP3 by gene
    z3e = entrez("ZP3")
    print("  ZP3 entrez:", z3e)
    # Pull ZP3 for all samples (chunked GET molecular-data, single gene)
    rows = {}
    chunks = [samples[i:i+150] for i in range(0,len(samples),150)]
    for ci,chk in enumerate(chunks):
        qs = "&".join(["sampleIds=%s"%s for s in chk])  # may be too long; use fetch POST
        # use POST fetch which supports bulk
        body={"molecularProfileIds":[mp],"sampleIds":chk,"entrezGeneIds":[z3e]}
        try:
            r=requests.post(API+"/molecular-data/fetch",json=body,headers=HEADERS,timeout=60)
            for rec in r.json():
                rows[rec["sampleId"]]=rec["value"]
        except Exception as e:
            print("  chunk %d err %s"%(ci,e))
        time.sleep(0.2)
    print("  Samples with ZP3 expression fetched:", len(rows))
    # 3) lineage/type metadata
    cts = get("/studies/ccle_broad_2019/samples", params={"pageSize":len(samples)})
    meta={}
    if isinstance(cts,list):
        for s in cts:
            meta[s["sampleId"]]=s
    # 4) Summary
    recrows=[]
    for sid,val in rows.items():
        recrows.append({"sampleId":sid,"ZP3":val,
                        "cancerType":meta.get(sid,{}).get("cancerTypeId",""),
                        "patientId":meta.get(sid,{}).get("patientId","")})
    df=pd.DataFrame(recrows)
    df.to_csv(os.path.join(BASE,"ccle_zp3_expr.csv"),index=False)
    print("  Saved ccle_zp3_expr.csv  shape:", df.shape)
    print("  Number of cell lines with ZP3>0:", int((df['ZP3']>0).sum()), "/", len(df))
    # 5) Key: CNS/glioma cell lines
    cns = df[df['cancerType'].astype(str).str.upper().str.contains("CNS",na=False)]
    print("  CNS lineage cell lines:",len(cns))
    if len(cns):
        print("  CNS with ZP3>0:",int((cns['ZP3']>0).sum()),"/",len(cns))
        print("  CNS ZP3 expression max/mean(>0):", cns['ZP3'].max(), cns.loc[cns['ZP3']>0,'ZP3'].mean() if (cns['ZP3']>0).any() else 0)
        print("  CNS ZP3 expression top8:")
        print(cns.sort_values('ZP3',ascending=False).head(8)[['sampleId','patientId','ZP3']].to_string(index=False))
    # 6) Comparison of other solid tumor cell lines (to see whether ZP3 is mainly in reproductive/other lineages)
    print("\n  Non-CNS ZP3 expression top5:")
    others = df[~df['cancerType'].astype(str).str.upper().str.contains("CNS",na=False)]
    print(others.sort_values('ZP3',ascending=False).head(5)[['sampleId','patientId','cancerType','ZP3']].to_string(index=False))
else:
    print("  CCLE sample-list fetch failed:", sl)

# ---------------- Test B: bulk pan-myeloid index ----------------
print("\n"+"="*70)
print("[Test B] bulk-level: whether ZP3↔TREM2 is independent of total myeloid burden")
marker_syms=["CD68","CD14","LYZ","CSF1R","ITGAM"]
for study in ["gbm_tcga","lgg_tcga"]:
    mp_r="%s_rna_seq_v2_mrna"%study
    sll="%s_rna_seq_v2_mrna"%study
    samples=get("/sample-lists/%s/sample-ids"%sll)
    if not isinstance(samples,list) or not samples:
        print("  %s sample list err"%study); continue
    # Merge into existing expr (already includes ZP3, TREM2, etc.)
    epath=os.path.join(BASE,"expr_%s_tcga_patient.csv"%("gbm" if "gbm" in study else "lgg"))
    expr=pd.read_csv(epath,index_col=0)
    # Need to align by patientId; expr index is patientId (previous script used patient level)
    # Map from sample to patient: fetch returns patientId
    for sym in marker_syms:
        e=entrez(sym)
        if not e: continue
        body={"molecularProfileIds":[mp_r],"sampleIds":samples[:600],"entrezGeneIds":[e]}
        try:
            r=requests.post(API+"/molecular-data/fetch",json=body,headers=HEADERS,timeout=60)
            for rec in r.json():
                pid=rec.get("patientId");val=rec.get("value")
                if pid and val is not None:
                    expr.loc[pid,sym]=float(val)
        except Exception as ex:
            print("  %s %s err %s"%(study,sym,ex))
        time.sleep(0.2)
    # Myeloid index = z-mean of available markers
    avail=[g for g in marker_syms if g in expr.columns and expr[g].notna().sum()>30]
    if avail:
        zs=(expr[avail]-expr[avail].mean())/expr[avail].std()
        mx=zs.mean(axis=1)
        expr["myeloid_idx"]=mx
        sub=expr.dropna(subset=["ZP3","TREM2","myeloid_idx"])
        import scipy.stats as st
        r0,p0=st.pearsonr(sub["ZP3"],sub["TREM2"])
        r1,p1=st.pearsonr(sub["ZP3"],sub["myeloid_idx"])
        r2,p2=st.pearsonr(sub["TREM2"],sub["myeloid_idx"])
        # partial corr
        def reg(y,x):
            X=np.column_stack([np.ones(len(x)),x]); c,_,_,_=np.linalg.lstsq(X,y,rcond=None)
            return c[0]+X[:,1]*c[1]
        rzp=sub["ZP3"].values-reg(sub["ZP3"].values,sub["myeloid_idx"].values)
        rtr=sub["TREM2"].values-reg(sub["TREM2"].values,sub["myeloid_idx"].values)
        rp,pp=st.pearsonr(rzp,rtr)
        print("--- %s n=%d ---"%(study,len(sub)))
        print("  ZP3 vs TREM2 crude: r=%.3f p=%.3g"%(r0,p0))
        print("  ZP3 vs myeloid index: r=%.3f p=%.3g"%(r1,p1))
        print("  TREM2 vs myeloid index: r=%.3f p=%.3g"%(r2,p2))
        print("  Partial correlation (controlling for myeloid index) ZP3 vs TREM2: r_partial=%.3f p=%.3g"%(rp,pp))
        expr.to_csv(os.path.join(BASE,"expr_%s_b_myeloid.csv"%("gbm" if "gbm" in study else "lgg")))
