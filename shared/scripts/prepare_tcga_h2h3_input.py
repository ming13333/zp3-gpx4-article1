# -*- coding: utf-8 -*-
"""
h2_h3_tcga.py adapter input generator
Background: h2_h3_tcga.py expects Xena-format input (HiSeq_TCGA_gene.xena.gz + GBM/LGG_clinicalMatrix.gz),
But the Xena S3 data source is currently unreachable (403). This script builds equivalent input from local real data:
  1. Expression matrix: stream-extract ZP3+ immune genes from 1.3GB TcgaTargetGtex_rsem_gene_tpm.gz,
     Filter TCGA-GBM/LGG tumor samples using GDC participant->cancer-type mapping (real data, not simulated)
  2. Clinical matrix: fetch OS_MONTHS/OS_STATUS for gbm_tcga/lgg_tcga from the cBioPortal API (real clinical data)
Outputs are written to output/h2_bulk/TCGA.GBM.sampleMap/ and TCGA.LGG.sampleMap/,
So h2_h3_tcga.py can be run end-to-end unchanged (without modifying the target script).

Usage: python prepare_tcga_h2h3_input.py
"""
import os, sys, gzip, json, time
import numpy as np
import pandas as pd
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
HUB = os.path.join(os.path.dirname(BASE), "phase1_knowledge_gap_filling")  # output/phase1_knowledge_gap_filling
TPM_GZ = os.path.join(HUB, "TcgaTargetGtex_rsem_gene_tpm.gz")
DISEASE_MAP = os.path.join(os.path.dirname(BASE), "tcga_pancan", "tcga_disease_map.json")
ENSG_MAP = os.path.join(os.path.dirname(BASE), "tcga_pancan", "ensg_map.json")

# Same gene set as h2_h3_tcga.py
IMMUNOSUPP_GENES = ["TGFB1", "IL10", "FOXP3", "CD274", "PDCD1", "CTLA4",
                    "MRC1", "CD163", "VSIG4", "ARG1", "IDO1", "VEGFA",
                    "CCL2", "CXCL12", "MSR1", "TREM2"]
M2_GENES = ["MRC1", "CD163", "MSR1", "ARG1", "TGFB1", "IL10", "VSIG4"]
TREG_GENES = ["FOXP3", "IL2RA", "CTLA4", "TIGIT"]
CHECKPT_GENES = ["CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3"]
GENES = sorted(set(["ZP3"] + IMMUNOSUPP_GENES + M2_GENES + TREG_GENES + CHECKPT_GENES))

CBIO = "https://www.cbioportal.org/api"

def http_get_json(url, timeout=60):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())

def build_expression():
    """Stream through TPM, extract target gene rows, filter samples by cancer type."""
    with open(DISEASE_MAP) as f:
        disease = json.load(f)          # submitter_id(TCGA-XX-XXXX) -> cancer
    with open(ENSG_MAP) as f:
        ensg_map = json.load(f)         # symbol -> ensg
    # Fill in target genes missing from ensg_map (official Ensembl GRCh38 IDs)
    EXTRA = {"VEGFA": "ENSG00000112715", "CCL2": "ENSG00000108691",
             "CXCL12": "ENSG00000107562", "TREM2": "ENSG00000095970"}
    ensg_map.update({k: v for k, v in EXTRA.items() if k in GENES})
    want = {v: k for k, v in ensg_map.items() if k in GENES}
    print(f"Target genes (Ensembl): {len(want)}")
    rows = {}                           # ensg -> {sample: value}
    with gzip.open(TPM_GZ, "rt") as f:
        header = f.readline().rstrip().split("\t")
        samples = header[1:]
        n = 0
        while True:
            lines = f.readlines(50000)
            if not lines:
                break
            for ln in lines:
                parts = ln.rstrip().split("\t")
                base = parts[0].split(".")[0]   # Remove Ensembl version suffix
                if base in want:
                    rows[base] = {s: float(v) for s, v in zip(samples, parts[1:])}
            n += len(lines)
            if n % 200000 == 0:
                print(f"  Scanned {n} lines, hit {len(rows)}")
    print(f"Scan complete, hit {len(rows)} gene rows")
    df = pd.DataFrame(rows)             # index=sample, columns=ensg
    # Filter TCGA tumor samples (01) + GBM/LGG
    df["cancer"] = df.index.map(lambda s: disease.get(s[:12], ""))
    gbm = df[(df["cancer"] == "GBM") & (df.index.str.contains("-01$"))]
    lgg = df[(df["cancer"] == "LGG") & (df.index.str.contains("-01$"))]
    print(f"  GBM tumor samples: {len(gbm)} | LGG tumor samples: {len(lgg)}")
    out = {}
    for cancer, sub in (("GBM", gbm), ("LGG", lgg)):
        d = sub.drop(columns=["cancer"]).T.rename(index=want)  # rows=symbol, cols=sample
        # Xena format: log2(TPM+1) is already the raw file value, just transpose and keep
        os.makedirs(os.path.join(BASE, f"TCGA.{cancer}.sampleMap"), exist_ok=True)
        p = os.path.join(BASE, f"TCGA.{cancer}.sampleMap", f"HiSeq_TCGA_gene.xena.gz")
        d.to_csv(p, sep="\t", compression="gzip")
        print(f"  Writing {p} ({d.shape[0]} genes x {d.shape[1]} samples)")
        out[cancer] = d
    return out

def build_clinical():
    """Pull OS_MONTHS/OS_STATUS from cBioPortal + sampleId->patientId mapping, construct clinicalMatrix."""
    for cancer in ("GBM", "LGG"):
        study = f"{cancer.lower()}_tcga"
        # sample -> patient
        samples = http_get_json(f"{CBIO}/studies/{study}/samples?projection=SUMMARY")
        sp = {s["sampleId"]: s["patientId"] for s in samples}
        # patient-level OS clinical
        clin = http_get_json(f"{CBIO}/studies/{study}/clinical-data?clinicalDataType=PATIENT&pageSize=10000")
        osd = {}
        for rec in clin:
            attr = rec.get("clinicalAttributeId")
            if attr in ("OS_MONTHS", "OS_STATUS"):
                osd.setdefault(rec["patientId"], {})[attr] = rec.get("value")
        # patient -> sample (take 01)
        rows = {}
        for pid, d in osd.items():
            sid = next((s for s, p in sp.items() if p == pid), None)
            if sid and "OS_MONTHS" in d and "OS_STATUS" in d:
                rows[sid] = {"OS_MONTHS": d["OS_MONTHS"], "OS_STATUS": d["OS_STATUS"]}
        df = pd.DataFrame.from_dict(rows, orient="index")
        p = os.path.join(BASE, f"TCGA.{cancer}.sampleMap", f"{cancer}_clinicalMatrix.gz")
        df.to_csv(p, sep="\t", compression="gzip")
        print(f"  Wrote {p} ({len(df)} samples)")

if __name__ == "__main__":
    print("=== 1/2 Building expression matrix (local real TPM) ===")
    build_expression()
    print("\n=== 2/2 Building clinical matrix (cBioPortal real OS) ===")
    build_clinical()
    print("\nInput ready. You can run h2_h3_tcga.py for end-to-end validation.")
