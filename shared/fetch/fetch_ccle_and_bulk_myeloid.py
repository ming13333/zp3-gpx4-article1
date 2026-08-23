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
方法3 + 检验B fetch：CCLE 纯肿瘤细胞系 ZP3 表达 + bulk 泛髓系指数偏相关
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

# ---------------- 方法3: CCLE ----------------
print("="*70)
print("[方法3] CCLE 纯肿瘤细胞系 ZP3 表达（检验:纯肿瘤细胞是否表达 ZP3）")
mp = "ccle_broad_2019_rna_seq_mrna"
# 1) sample list（全部 CCLE cell lines）
sl = get("/sample-lists/ccle_broad_2019_all/sample-ids")
if isinstance(sl, list) and sl:
    samples = sl
    print("  CCLE 全部细胞系:", len(samples))
    # 2) 逐基因拉 ZP3
    z3e = entrez("ZP3")
    print("  ZP3 entrez:", z3e)
    # 拉 ZP3 全部样本(分块 GET molecular-data, 单基因)
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
    print("  拉到 ZP3 表达样本数:", len(rows))
    # 3) lineage/type 元数据
    cts = get("/studies/ccle_broad_2019/samples", params={"pageSize":len(samples)})
    meta={}
    if isinstance(cts,list):
        for s in cts:
            meta[s["sampleId"]]=s
    # 4) 汇总
    recrows=[]
    for sid,val in rows.items():
        recrows.append({"sampleId":sid,"ZP3":val,
                        "cancerType":meta.get(sid,{}).get("cancerTypeId",""),
                        "patientId":meta.get(sid,{}).get("patientId","")})
    df=pd.DataFrame(recrows)
    df.to_csv(os.path.join(BASE,"ccle_zp3_expr.csv"),index=False)
    print("  已存 ccle_zp3_expr.csv  形状:",df.shape)
    print("  ZP3>0 的细胞系数:",int((df['ZP3']>0).sum()),"/",len(df))
    # 5) 重点: CNS/胶质瘤细胞系
    cns = df[df['cancerType'].astype(str).str.upper().str.contains("CNS",na=False)]
    print("  CNS lineage 细胞系:",len(cns))
    if len(cns):
        print("  CNS 中 ZP3>0:",int((cns['ZP3']>0).sum()),"/",len(cns))
        print("  CNS ZP3 表达 max/mean(>0):", cns['ZP3'].max(), cns.loc[cns['ZP3']>0,'ZP3'].mean() if (cns['ZP3']>0).any() else 0)
        print("  CNS ZP3 表达 top8:")
        print(cns.sort_values('ZP3',ascending=False).head(8)[['sampleId','patientId','ZP3']].to_string(index=False))
    # 6) 其它实体瘤细胞系对比（看 ZP3 是否主要在生殖/其它谱系）
    print("\n  非中枢 ZP3 表达 top5:")
    others = df[~df['cancerType'].astype(str).str.upper().str.contains("CNS",na=False)]
    print(others.sort_values('ZP3',ascending=False).head(5)[['sampleId','patientId','cancerType','ZP3']].to_string(index=False))
else:
    print("  CCLE sample-list 拉取失败:", sl)

# ---------------- 检验B: bulk 泛髓系指数 ----------------
print("\n"+"="*70)
print("[检验 B] bulk 层面 ZP3↔TREM2 是否独立于总髓系负荷")
marker_syms=["CD68","CD14","LYZ","CSF1R","ITGAM"]
for study in ["gbm_tcga","lgg_tcga"]:
    mp_r="%s_rna_seq_v2_mrna"%study
    sll="%s_rna_seq_v2_mrna"%study
    samples=get("/sample-lists/%s/sample-ids"%sll)
    if not isinstance(samples,list) or not samples:
        print("  %s sample list err"%study); continue
    # 并入现有 expr（已含 ZP3,TREM2 等）
    epath=os.path.join(BASE,"expr_%s_tcga_patient.csv"%("gbm" if "gbm" in study else "lgg"))
    expr=pd.read_csv(epath,index_col=0)
    # 需要按 patientId 对齐；expr index 就是 patientId（之前脚本用 patient 级）
    # 从 sample 拉到 patient: fetch 返回含 patientId
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
    # 髓系指数 = 可用 marker 的 z-mean
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
        print("  ZP3vsTREM2 粗: r=%.3f p=%.3g"%(r0,p0))
        print("  ZP3vs髓系指数: r=%.3f p=%.3g"%(r1,p1))
        print("  TREM2vs髓系指数: r=%.3f p=%.3g"%(r2,p2))
        print("  偏相关(控髓系指数) ZP3vsTREM2: r_partial=%.3f p=%.3g"%(rp,pp))
        expr.to_csv(os.path.join(BASE,"expr_%s_b_myeloid.csv"%("gbm" if "gbm" in study else "lgg")))
