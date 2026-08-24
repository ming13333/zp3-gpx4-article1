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
Method 3: ZP3 expression in CCLE pure tumor cell lines (GET endpoint version)
Lineage: inferred from sampleId suffix (e.g., _CENTRAL_NERVOUS_SYSTEM)
"""
import requests, os, time
import numpy as np, pandas as pd
API="https://www.cbioportal.org/api"
BASE=os.path.join(ROOT, "output", "h2_bulk")
H={"Accept":"application/json"}
mp="ccle_broad_2019_rna_seq_mrna"

def get(url, params=None, retries=3):
    for i in range(retries):
        try:
            r=requests.get(url,params=params,headers=H,timeout=120)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if i==retries-1:
                print("ERR",e); return []
            time.sleep(2)

# Fetch full ZP3 data (GET, single request)
zs=get("%s/molecular-profiles/%s/molecular-data"%(API,mp),
       params={"sampleListId":"ccle_broad_2019_all","entrezGeneId":7784})
print("ZP3 records fetched:", len(zs))
rows=[]
for r in zs:
    rows.append({"sampleId":r.get("sampleId"),"ZP3":r.get("value")})
df=pd.DataFrame(rows)
df.to_csv(os.path.join(BASE,"ccle_zp3_expr.csv"),index=False)
print("Saved ccle_zp3_expr.csv", df.shape)

# lineage from sampleId suffix
def lineage(sid):
    sid=str(sid)
    if "_" not in sid: return "other"
    return sid.split("_",1)[1] if sid.count("_")>=1 else "other"
df["lineage"]=df["sampleId"].apply(lineage)

print("\n=== Overview ===")
print("Total cell lines:", len(df), "| ZP3>0:", int((df['ZP3']>0).sum()), "| expression rate=%.1f%%"%(100*(df['ZP3']>0).mean()))
print("\n=== CNS lineage (glioma and other CNS) cell lines ZP3 ===")
cns=df[df['lineage'].astype(str).str.upper()=="CENTRAL_NERVOUS_SYSTEM"]
print("CNS cell line count:", len(cns), "| ZP3>0:", int((cns['ZP3']>0).sum()), "| expression rate=%.1f%%"%(100*(cns['ZP3']>0).mean()))
print("CNS ZP3 expression max=%.2f, mean(all)=%.3f, mean(>0)=%.3f"%(cns['ZP3'].max(), cns['ZP3'].mean(), cns.loc[cns['ZP3']>0,'ZP3'].mean() if (cns['ZP3']>0).any() else 0))
print("CNS expression top10:")
print(cns.sort_values('ZP3',ascending=False).head(10)[['sampleId','ZP3']].to_string(index=False))

print("\n=== Non-CNS lineage background ===")
other=df[df['lineage'].astype(str).str.upper()!="CENTRAL_NERVOUS_SYSTEM"]
print("Non-CNS cell lines:", len(other), "| ZP3>0:", int((other['ZP3']>0).sum()), "| expr rate=%.1f%%"%(100*(other['ZP3']>0).mean()))
print("Non-CNS expression rate top lineages:")
top=df.groupby('lineage')['ZP3'].agg(['count',lambda s:(s>0).mean(),'mean','max']).rename(columns={'count':'n','<lambda_0>':'frac_pos','mean':'mean','max':'max'})
print(top.sort_values('mean',ascending=False).head(12).to_string())

print("\n=== Reproductive/ovary and other classic ZP3 lineages (controls) ===")
for lin in ['OVARY','TESTIS','ENDOMETRIUM','CERVIX','BREAST']:
    s=df[df['lineage'].astype(str).str.upper()==lin]
    if len(s):
        print("  %-12s n=%3d ZP3>0=%3d expr rate=%5.1f%% mean=%.2f"%(lin,len(s),int((s['ZP3']>0).sum()),100*(s['ZP3']>0).mean(),s['ZP3'].mean()))
