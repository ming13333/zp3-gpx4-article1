# -*- coding: utf-8 -*-
import requests, time
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
H={"User-Agent":"workbuddy/1.0"}
def fetch(pmids):
    p={"db":"pubmed","id":",".join(pmids),"rettype":"abstract","retmode":"text"}
    r=requests.get(BASE+"efetch.fcgi",params=p,headers=H,timeout=40)
    return r.text
for pid in ["39039845","38125942","40684330","41494530"]:
    print("\n"+"="*80)
    print("PMID:",pid)
    try:
        t=fetch([pid])
        print(t[:3200])
    except Exception as e:
        print("ERROR",e)
    time.sleep(0.5)
