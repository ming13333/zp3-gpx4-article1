# -*- coding: utf-8 -*-
"""PubMed E-utilities key verification search"""
import requests, time, json, sys

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
HEADERS = {"User-Agent": "workbuddy-research/1.0"}

def esearch(term, retmax=15):
    p = {"db":"pubmed","term":term,"retmax":retmax,"retmode":"json"}
    r = requests.get(BASE+"esearch.fcgi", params=p, headers=HEADERS, timeout=30)
    return r.json()["esearchresult"]

def efetch(pmids):
    if not pmids: return []
    p = {"db":"pubmed","id":",".join(pmids),"rettype":"medline","retmode":"text"}
    r = requests.get(BASE+"efetch.fcgi", params=p, headers=HEADERS, timeout=40)
    return [b for b in r.text.split("\n\n") if b.strip()]

def parse_medline(block):
    d={}
    for line in block.splitlines():
        if line.startswith("TI  -"): d.setdefault("title",[]).append(line[6:].strip())
        elif line.startswith("JT  -"): d["journal"]=line[6:].strip()
        elif line.startswith("DP  -"): d["date"]=line[6:].strip()
        elif line.startswith("PMID"): d["pmid"]=line[6:].strip()
        elif line.startswith("AB  -"): d.setdefault("ab",[]).append(line[6:].strip())
    return d

QUERIES = {
 "Q1_ZP3_tumor_DC": 'ZP3[Title/Abstract] AND (dendritic OR DC OR "tumor microenvironment")',
 "Q2_ZP3_CNS": 'ZP3[Title/Abstract] AND (brain OR glioma OR microglia OR neuro OR astrocytoma)',
 "Q3_ZP3_singlecell_immune": 'ZP3 AND (single-cell OR scRNA-seq OR myeloid OR macrophage) AND 2019:2026[dp]',
 "Q4_extracellular_GPX4_neuro": 'GPX4[Title/Abstract] AND (brain OR stroke OR glioma OR ischemia OR neurological) AND 2024:2026[dp]',
 "Q5_ZP3_Cancer_transcript": '(ZP3-Cancer OR "ZP3-Cancer" OR "Zona Pellucida glycoprotein 3" OR ZP3) AND (cancer antigen OR tumor antigen OR ectopic expression) AND 2020:2026[dp]',
}

out=[]
for name,q in QUERIES.items():
    try:
        res=esearch(q)
        n=int(res.get("count",0)); ids=res.get("idlist",[])
        out.append(f"\n===== {name} | count={n} =====")
        if ids:
            recs=efetch(ids[:8])
            for b in recs:
                m=parse_medline(b)
                out.append(f"[{m.get('pmid','')}] {m.get('title',[''])[0]} | {m.get('journal','')} | {m.get('date','')}")
    except Exception as e:
        out.append(f"{name} ERROR: {e}")
    time.sleep(0.4)

print("\n".join(out))
