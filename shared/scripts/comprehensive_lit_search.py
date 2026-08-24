# -*- coding: utf-8 -*-
"""
Cross-database literature search: PubMed (E-utilities) + Europe PMC REST
Output RIS format, including PMID/DOI/Title/Abstract/Journal/Year/Keywords
"""
import os, requests, time, json, re
from datetime import datetime

def _project_root():
    d = os.path.dirname(os.path.abspath(__file__))
    while True:
        if os.path.isdir(os.path.join(d, "output")):
            return d
        p = os.path.dirname(d)
        if p == d:
            break
        d = p
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = _project_root()

BASE_PUB = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BASE_EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"
HEADERS = {"User-Agent": "workbuddy-research/1.0"}

# Core search strategy (grouped by topic)
QUERIES = {
    # === Core concept: Cell 2026 GPX4-ZP3 axis ===
    "core_GPX4_ZP3_Cell2026": '("Extracellular GPX4" OR "GPX4" AND "ZP3") AND (dendritic OR immunity OR antitumor) AND 2025:2026[dp]',
    
    # === Ectopic expression and function of ZP3 in tumors/cancer ===
    "ZP3_cancer_ectopic": '(ZP3[Title/Abstract] OR "Zona Pellucida Glycoprotein 3"[Title/Abstract]) AND (cancer OR carcinoma OR tumor OR neoplasm) AND 2018:2026[dp]',
    
    # === ZP3-Cancer alternative transcript ===
    "ZP3_Cancer_transcript": '("ZP3-Cancer" OR "ZP3 Cancer" OR "alternative transcript" ZP3) AND 2020:2026[dp]',
    
    # === ZP3 and Notch/migration/adhesion signaling ===
    "ZP3_Notch_signaling": '(ZP3 AND (Notch OR migration OR adhesion OR invasion)) AND 2018:2026[dp]',
    
    # === Extracellular GPX4 / ferroptosis DAMP / immunity ===
    "extracellular_GPX4_DAMP": '("extracellular GPX4" OR "GPX4 DAMP" OR "ferroptosis DAMP" OR "GPX4 release") AND 2020:2026[dp]',
    
    # === Neuro field: glioma/brain tumor single-cell myeloid/immunity ===
    "glioma_scRNA_myeloid": '(glioma OR glioblastoma OR GBM) AND (single-cell OR scRNA-seq OR "single nucleus") AND (myeloid OR macrophage OR microglia OR "tumor-associated") AND 2020:2026[dp]',
    
    # === Stroke/intracerebral hemorrhage ferroptosis immunity ===
    "stroke_ICH_ferroptosis_immunity": '(stroke OR "intracerebral hemorrhage" OR ICH OR "cerebral ischemia") AND (ferroptosis OR GPX4) AND (immune OR inflammation OR microglia OR macrophage) AND 2020:2026[dp]',
    
    # === ZP3 expression in neural/brain tissue (eQTL/biomarkers) ===
    "ZP3_brain_expression": '(ZP3 AND (brain OR cerebral OR neurological OR "PsychENCODE" OR "GTEx" OR eQTL)) AND 2018:2026[dp]',
    
    # === Ferroptosis in neurodegeneration/neuroinflammation ===
    "ferroptosis_neurodegeneration": '(ferroptosis AND (Parkinson OR Alzheimer OR ALS OR "multiple sclerosis" OR neuroinflammation OR neurodegeneration)) AND 2020:2026[dp]',
    
    # === DC/myeloid cell cAMP-PRKA glycolysis immunosuppression ===
    "DC_cAMP_glycolysis_immunity": '(dendritic AND (cAMP OR PRKA OR PKA) AND glycolysis AND (maturation OR immunosuppression)) AND 2018:2026[dp]',
}

def pubmed_esearch(term, retmax=50):
    p = {"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json", "usehistory": "y"}
    r = requests.get(BASE_PUB + "esearch.fcgi", params=p, headers=HEADERS, timeout=40)
    return r.json()["esearchresult"]

def pubmed_efetch(pmids, batch=200):
    if not pmids:
        return []
    out = []
    for i in range(0, len(pmids), batch):
        batch_ids = pmids[i:i+batch]
        p = {"db": "pubmed", "id": ",".join(batch_ids), "rettype": "medline", "retmode": "text"}
        r = requests.get(BASE_PUB + "efetch.fcgi", params=p, headers=HEADERS, timeout=60)
        out.extend([b for b in r.text.split("\n\n") if b.strip()])
        time.sleep(0.35)
    return out

def parse_medline(block):
    d = {"pmid": "", "doi": "", "title": "", "journal": "", "year": "", "abstract": "", "authors": [], "keywords": []}
    lines = block.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("PMID- "):
            d["pmid"] = line[6:].strip()
        elif line.startswith("TI  - "):
            d["title"] = line[6:].strip()
            i += 1
            while i < len(lines) and lines[i].startswith("      "):
                d["title"] += " " + lines[i].strip()
                i += 1
            i -= 1
        elif line.startswith("JT  - "):
            d["journal"] = line[6:].strip()
        elif line.startswith("DP  - "):
            d["year"] = line[6:].strip()[:4]
        elif line.startswith("AB  - "):
            d["abstract"] = line[6:].strip()
            i += 1
            while i < len(lines) and lines[i].startswith("      "):
                d["abstract"] += " " + lines[i].strip()
                i += 1
            i -= 1
        elif line.startswith("FAU - "):
            d["authors"].append(line[6:].strip())
        elif line.startswith("AID - ") and "[doi]" in line:
            d["doi"] = line[6:].split("[doi]")[0].strip()
        elif line.startswith("OT  - "):
            d["keywords"].append(line[6:].strip())
        elif line.startswith("MH  - "):
            d["keywords"].append(line[6:].strip())
        i += 1
    return d

def europepmc_search(query, page_size=50):
    """Search Europe PMC REST API, return JSON"""
    url = BASE_EPMC + "search"
    params = {"query": query, "format": "json", "pageSize": page_size, "resultType": "lite"}
    try:
        r = requests.get(url, params=params, headers=HEADERS, timeout=40)
        return r.json()
    except:
        return {}

def europepmc_fetch(pmids):
    """Batch fetch Europe PMC details"""
    if not pmids:
        return []
    url = BASE_EPMC + "search"
    id_query = " OR ".join([f"EXT_ID:{p}" for p in pmids])
    params = {"query": id_query, "format": "json", "pageSize": len(pmids), "resultType": "lite"}
    try:
        r = requests.get(url, params=params, headers=HEADERS, timeout=40)
        return r.json().get("resultList", {}).get("result", [])
    except:
        return []

def parse_europepmc(rec):
    d = {"pmid": rec.get("pmid", ""), "doi": rec.get("doi", ""), "title": rec.get("title", ""),
         "journal": rec.get("journalInfo", {}).get("journal", {}).get("title", ""),
         "year": str(rec.get("pubYear", "")), "abstract": rec.get("abstractText", ""),
         "authors": [a.get("fullName", "") for a in rec.get("authorList", {}).get("author", [])],
         "keywords": []}
    return d

def to_ris(record):
    """Convert a dict to a RIS entry"""
    lines = ["TY  - JOUR"]
    if record.get("pmid"): lines.append(f"ID  - PMID:{record['pmid']}")
    if record.get("doi"): lines.append(f"DO  - {record['doi']}")
    if record.get("title"): lines.append(f"TI  - {record['title']}")
    for au in record.get("authors", []):
        lines.append(f"AU  - {au}")
    if record.get("journal"): lines.append(f"JO  - {record['journal']}")
    if record.get("year"): lines.append(f"PY  - {record['year']}")
    if record.get("abstract"): lines.append(f"AB  - {record['abstract']}")
    for kw in record.get("keywords", []):
        lines.append(f"KW  - {kw}")
    lines.append("ER  - ")
    return "\n".join(lines)

def main():
    all_records = []
    seen_pmids = set()
    seen_dois = set()
    
    for name, query in QUERIES.items():
        print(f"\n=== {name} ===")
        
        # 1) PubMed
        try:
            res = pubmed_esearch(query)
            count = int(res.get("count", 0))
            ids = res.get("idlist", [])
            print(f"  PubMed: {count} hits, fetching top {len(ids)}")
            if ids:
                blocks = pubmed_efetch(ids)
                for b in blocks:
                    rec = parse_medline(b)
                    if rec["pmid"] and rec["pmid"] not in seen_pmids:
                        seen_pmids.add(rec["pmid"])
                        all_records.append(rec)
                        if rec.get("doi"): seen_dois.add(rec["doi"])
        except Exception as e:
            print(f"  PubMed error: {e}")
        
        # 2) Europe PMC (supplement & deduplicate)
        try:
            epmc = europepmc_search(query)
            hits = epmc.get("hitCount", 0)
            results = epmc.get("resultList", {}).get("result", [])
            print(f"  Europe PMC: {hits} hits, got {len(results)} records")
            for rec in results:
                p = parse_europepmc(rec)
                pmid = p.get("pmid", "")
                doi = p.get("doi", "")
                if (pmid and pmid not in seen_pmids) or (doi and doi not in seen_dois):
                    if pmid: seen_pmids.add(pmid)
                    if doi: seen_dois.add(doi)
                    all_records.append(p)
        except Exception as e:
            print(f"  Europe PMC error: {e}")
        
        time.sleep(0.5)
    
    print(f"\n=== Total after deduplication: {len(all_records)} papers ===")
    
    # Write RIS
    ris_path = os.path.join(ROOT, "output", "literature_background.ris")
    with open(ris_path, "w", encoding="utf-8") as f:
        for rec in all_records:
            f.write(to_ris(rec) + "\n\n")
    print(f"RIS written: {ris_path}")
    
    # Also output CSV for easy browsing
    import csv
    csv_path = os.path.join(ROOT, "output", "literature_background.csv")
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["PMID", "DOI", "Title", "Journal", "Year", "FirstAuthor", "Keywords", "Abstract"])
        for rec in all_records:
            w.writerow([
                rec.get("pmid", ""), rec.get("doi", ""), rec.get("title", ""),
                rec.get("journal", ""), rec.get("year", ""),
                rec.get("authors", [""])[0] if rec.get("authors") else "",
                "; ".join(rec.get("keywords", [])),
                rec.get("abstract", "")[:500]
            ])
    print(f"CSV written: {csv_path}")

if __name__ == "__main__":
    main()
