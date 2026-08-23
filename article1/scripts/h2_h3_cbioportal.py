# -*- coding: utf-8 -*-
"""
H2 (预后关联) + H3 (免疫抑制 TME 关联) — TCGA GBM / LGG 泛胶质瘤 bulk RNA-seq
数据来源：cBioPortal REST API v2 (https://www.cbioportal.org/api)
  - clinical-data (PATIENT): OS_MONTHS / OS_STATUS
  - molecular-data/fetch: RNA-seq V2 RSEM (gbm_tcga_rna_seq_v2_mrna, lgg_tcga_rna_seq_v2_mrna)
描述级分析，诚实标注：bulk RNA 无法区分 ZP3 经典转录本 vs ZP3-Cancer 替代异构体。
"""
import os, time, sys
import numpy as np
import pandas as pd
import requests
from scipy import stats

API = "https://www.cbioportal.org/api"
BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output", "h2_bulk")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article1", "results")

# 统一使用经 lifelines/标准实现验证的 log-rank（消除旧版两套不一致手写实现）
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "shared", "scripts"))
from stats_utils import logrank  # noqa: E402

# 免疫抑制相关标志基因（文献常见，作描述性关联）
IMMUNOSUPP_GENES = ["TGFB1","IL10","FOXP3","CD274","PDCD1","CTLA4","MRC1","CD163",
                    "VSIG4","ARG1","IDO1","VEGFA","CCL2","CXCL12","MSR1","TREM2"]
M2_GENES    = ["MRC1","CD163","MSR1","ARG1","TGFB1","IL10","VSIG4"]
TREG_GENES  = ["FOXP3","IL2RA","CTLA4","TIGIT"]
CHECKPT_GENES = ["CD274","PDCD1","CTLA4","HAVCR2","LAG3"]

SYM2ENT = {}  # 全局：symbol -> entrez

def api_get(path, params=None):
    r = requests.get(API + path, params=params or {}, timeout=60)
    r.raise_for_status()
    return r.json()

def api_get_all(path, params=None):
    """带分页的列表拉取（cBioPortal 通过 nextPageToken 响应头翻页）"""
    params = dict(params or {})
    out = []
    token = None
    while True:
        p = dict(params)
        if token:
            p["pageToken"] = token
        r = requests.get(API + path, params=p, timeout=60)
        r.raise_for_status()
        d = r.json()
        if not isinstance(d, list):
            return d
        out.extend(d)
        nt = r.headers.get("nextPageToken") or r.headers.get("X-Next-Page-Token")
        if not nt:
            break
        token = nt
    return out

def api_post(path, body):
    r = requests.post(API + path, json=body, timeout=120)
    r.raise_for_status()
    return r.json()

def resolve_entrez(symbols):
    m = {}
    for s in symbols:
        try:
            d = api_get(f"/genes/{s}")
            m[s] = d.get("entrezGeneId")
        except Exception as e:
            print(f"  !! 无法解析 {s}: {e}")
    return m

def get_sample_ids(study):
    return api_get(f"/sample-lists/{study}_all/sample-ids")

def get_clinical_os(study):
    """返回 patientId -> (os_time_months, os_event 0/1) 仅含 OS_MONTHS/OS_STATUS"""
    data = api_get_all(f"/studies/{study}/clinical-data",
                       params={"clinicalDataType": "PATIENT", "pageSize": 100000})
    pat = {}
    for rec in data:
        pid = rec.get("patientId"); attr = rec.get("clinicalAttributeId"); val = rec.get("value")
        if pid is None or attr is None:
            continue
        pat.setdefault(pid, {})[attr] = val
    os_time, os_event = {}, {}
    for pid, attrs in pat.items():
        t = attrs.get("OS_MONTHS")
        s = str(attrs.get("OS_STATUS", "")).upper()
        if t is None:
            continue
        try:
            t = float(t)
        except Exception:
            continue
        if t <= 0:
            continue
        if "DECEASE" in s or s.strip() == "1":
            ev = 1
        elif "LIV" in s or s.strip() == "0":
            ev = 0
        else:
            ev = np.nan
        os_time[pid] = t
        os_event[pid] = ev
    return os_time, os_event

def get_expression(study, profile, sample_list_id, entrez_list):
    """patient 级基因矩阵。优先【批量 POST】/molecular-profiles/{profile}/molecular-data/fetch
    （I1 修复：单次请求拉取全部 entrez，避免逐基因 GET；实测公共实例要求 profile 在 URL 路径、
    body 为单对象 {entrezGeneIds, sampleListId}，旧式 /molecular-data/fetch+数组 返回空），
    失败则回退逐基因。"""
    ent2sym = {v: k for k, v in SYM2ENT.items()}
    body = {
        "sampleListId": sample_list_id,
        "entrezGeneIds": [int(e) for e in entrez_list],
    }
    try:
        data = api_post(f"/molecular-profiles/{profile}/molecular-data/fetch", body)
        rows = {}
        for rec in data:
            pid = rec.get("patientId"); val = rec.get("value")
            if pid is None or val is None:
                continue
            sym = ent2sym.get(rec.get("entrezGeneId"))
            if sym is None:
                continue
            rows.setdefault(pid, {})[sym] = float(val)
        if rows:
            return pd.DataFrame.from_dict(rows, orient="index")
        print("    !! 批量返回为空，回退逐基因")
    except Exception as e:
        print(f"    !! 批量拉取失败，回退逐基因: {e}")
    # 回退：逐基因 GET
    rows = {}
    for ent in entrez_list:
        try:
            data = api_get(f"/molecular-profiles/{profile}/molecular-data",
                           params={"sampleListId": sample_list_id, "entrezGeneId": ent})
        except Exception as e:
            print(f"    !! 基因 entrez={ent} 拉取失败: {e}")
            continue
        for rec in data:
            pid = rec.get("patientId"); val = rec.get("value")
            if pid is None or val is None:
                continue
            sym = ent2sym.get(ent)
            if sym is None:
                continue
            rows.setdefault(pid, {})[sym] = float(val)
        time.sleep(0.05)
    return pd.DataFrame.from_dict(rows, orient="index")

def analyze(study, profile, genes_syms):
    print("\n" + "=" * 72)
    print(f"队列: {study} (profile={profile})")
    sample_list_id = f"{study}_rna_seq_v2_mrna"
    os_time, os_event = get_clinical_os(study)
    print(f"  有 OS 的 patient 数: {len(os_time)}")
    entrez = [SYM2ENT[s] for s in genes_syms if SYM2ENT.get(s)]
    expr = get_expression(study, profile, sample_list_id, entrez)
    print(f"  表达矩阵: {expr.shape[0]} patient × {expr.shape[1]} gene")
    if "ZP3" not in expr.columns:
        print("  !! 无 ZP3, 跳过"); return None
    zp3 = expr["ZP3"]
    time_s = pd.Series(os_time); event_s = pd.Series(os_event)
    merged = pd.DataFrame({"ZP3": zp3, "time": time_s, "event": event_s}).dropna()
    merged = merged[merged["time"] > 0]
    print(f"  H2 可用样本(表达+OS): {len(merged)}")
    if len(merged) < 30:
        print("  !! 样本不足"); return None
    med = merged["ZP3"].median()
    merged["group"] = (merged["ZP3"] > med).astype(int)
    hi = merged[merged.group == 1]; lo = merged[merged.group == 0]
    chi2, p = logrank(merged["time"].values, merged["event"].values, merged["group"].values)
    rate_hi = hi["event"].mean(); rate_lo = lo["event"].mean()
    direction = "High ZP3 预后更差" if rate_hi > rate_lo else "High ZP3 预后更好"
    print(f"  ZP3 中位={med:.3f} | High n={len(hi)} vs Low n={len(lo)}")
    print(f"  H2 logrank: chi2={chi2:.3f}, p={p:.4g}")
    print(f"  事件率 High={rate_hi:.3f} Low={rate_lo:.3f} -> {direction}")
    # H3: ZP3 vs 免疫抑制标志
    # 注：表达量为右侧偏态、非正态，改用 Spearman 秩相关（原为 Pearson，已修正）
    print("  --- H3: ZP3 vs 免疫抑制标志 (Spearman rho) ---")
    rows = []
    h3_genes = list(dict.fromkeys(IMMUNOSUPP_GENES + M2_GENES + TREG_GENES + CHECKPT_GENES))
    for gene in h3_genes:
        if gene in expr.columns:
            sub = pd.concat([merged["ZP3"], expr[gene]], axis=1).dropna()
            if len(sub) > 20:
                r, pp = stats.spearmanr(sub["ZP3"], sub[gene])
                rows.append((gene, round(float(r), 3), round(float(pp), 4), len(sub)))
    h3 = pd.DataFrame(rows, columns=["gene", "spearman_rho", "p", "n"]).sort_values("spearman_rho", ascending=False)
    print(h3.to_string(index=False))
    # 复合免疫抑制评分（M2+TREG+CHECKPT 可用基因 z-mean）
    sig = [g for g in M2_GENES + TREG_GENES + CHECKPT_GENES if g in expr.columns]
    immuno_r = np.nan; immuno_p = np.nan; immuno_n = 0
    if sig:
        sub = expr[sig].loc[merged.index]
        z = (sub - sub.mean()) / sub.std()
        immuno = z.mean(axis=1)
        gg = pd.concat([merged["ZP3"], immuno.rename("immuno_score")], axis=1).dropna()
        immuno_r, immuno_p = stats.spearmanr(gg["ZP3"], gg["immuno_score"])
        immuno_n = len(gg)
        print(f"  复合免疫抑制评分 vs ZP3: rho={immuno_r:.3f}, p={immuno_p:.4g}, n={immuno_n}")
    # 保存
    merged.to_csv(os.path.join(OUT, f"h2_{study}_zp3_os.csv"))
    h3.to_csv(os.path.join(OUT, f"h3_{study}_zp3_immuno.csv"), index=False)
    expr.to_csv(os.path.join(OUT, f"expr_{study}_patient.csv"))
    return {"study": study, "n": len(merged), "med": float(med), "chi2": float(chi2),
            "p": float(p), "dir": direction, "immuno_r": float(immuno_r) if sig else None,
            "immuno_p": float(immuno_p) if sig else None, "h3": h3}

if __name__ == "__main__":
    all_syms = list(dict.fromkeys(["ZP3"] + IMMUNOSUPP_GENES + M2_GENES + TREG_GENES + CHECKPT_GENES))
    print("解析基因 Entrez ID ...")
    SYM2ENT = resolve_entrez(all_syms)
    print("  ", {k: v for k, v in SYM2ENT.items()})
    results = {}
    results["gbm"] = analyze("gbm_tcga", "gbm_tcga_rna_seq_v2_mrna", all_syms)
    results["lgg"] = analyze("lgg_tcga", "lgg_tcga_rna_seq_v2_mrna", all_syms)
    print("\n=== H2/H3 分析完成（描述级，需外部验证；ZP3 异构体未区分）===")
