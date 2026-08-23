# -*- coding: utf-8 -*-
"""
h2_h3_tcga.py 的适配输入生成器
背景：h2_h3_tcga.py 期望 Xena 格式输入（HiSeq_TCGA_gene.xena.gz + GBM/LGG_clinicalMatrix.gz），
但 Xena S3 数据源当前 403 不可达。本脚本用【本地真实数据】构造等价输入：
  1. 表达矩阵：从 1.3GB TcgaTargetGtex_rsem_gene_tpm.gz 流式提取 ZP3+免疫基因，
     用 GDC 参与者→癌种映射筛选 TCGA-GBM/LGG 肿瘤样本（真实数据，非模拟）
  2. 临床矩阵：从 cBioPortal API 拉取 gbm_tcga/lgg_tcga 的 OS_MONTHS/OS_STATUS（真实临床数据）
产物写入 output/h2_bulk/TCGA.GBM.sampleMap/ 与 TCGA.LGG.sampleMap/，
使 h2_h3_tcga.py 可原样端到端运行（不改动目标脚本）。

用法: python prepare_tcga_h2h3_input.py
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

# 与 h2_h3_tcga.py 相同的基因集合
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
    """流式扫描 TPM，提取目标基因行，按癌种筛样本。"""
    with open(DISEASE_MAP) as f:
        disease = json.load(f)          # submitter_id(TCGA-XX-XXXX) -> cancer
    with open(ENSG_MAP) as f:
        ensg_map = json.load(f)         # symbol -> ensg
    # 补全 ensg_map 中缺失的目标基因（Ensembl GRCh38 官方 ID）
    EXTRA = {"VEGFA": "ENSG00000112715", "CCL2": "ENSG00000108691",
             "CXCL12": "ENSG00000107562", "TREM2": "ENSG00000095970"}
    ensg_map.update({k: v for k, v in EXTRA.items() if k in GENES})
    want = {v: k for k, v in ensg_map.items() if k in GENES}
    print(f"目标基因(Ensembl): {len(want)} 个")
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
                base = parts[0].split(".")[0]   # 去掉 Ensembl 版本号后缀
                if base in want:
                    rows[base] = {s: float(v) for s, v in zip(samples, parts[1:])}
            n += len(lines)
            if n % 200000 == 0:
                print(f"  已扫描 {n} 行, 命中 {len(rows)}")
    print(f"扫描完成，命中 {len(rows)} 个基因行")
    df = pd.DataFrame(rows)             # index=sample, columns=ensg
    # 筛 TCGA 肿瘤样本(01) + GBM/LGG
    df["cancer"] = df.index.map(lambda s: disease.get(s[:12], ""))
    gbm = df[(df["cancer"] == "GBM") & (df.index.str.contains("-01$"))]
    lgg = df[(df["cancer"] == "LGG") & (df.index.str.contains("-01$"))]
    print(f"  GBM 肿瘤样本: {len(gbm)} | LGG 肿瘤样本: {len(lgg)}")
    out = {}
    for cancer, sub in (("GBM", gbm), ("LGG", lgg)):
        d = sub.drop(columns=["cancer"]).T.rename(index=want)  # rows=symbol, cols=sample
        # Xena 格式：log2(TPM+1) 已是文件原始值，直接转置保留
        os.makedirs(os.path.join(BASE, f"TCGA.{cancer}.sampleMap"), exist_ok=True)
        p = os.path.join(BASE, f"TCGA.{cancer}.sampleMap", f"HiSeq_TCGA_gene.xena.gz")
        d.to_csv(p, sep="\t", compression="gzip")
        print(f"  写出 {p} ({d.shape[0]} 基因 x {d.shape[1]} 样本)")
        out[cancer] = d
    return out

def build_clinical():
    """从 cBioPortal 拉 OS_MONTHS/OS_STATUS + sampleId->patientId 映射，构造 clinicalMatrix。"""
    for cancer in ("GBM", "LGG"):
        study = f"{cancer.lower()}_tcga"
        # 样本 -> patient
        samples = http_get_json(f"{CBIO}/studies/{study}/samples?projection=SUMMARY")
        sp = {s["sampleId"]: s["patientId"] for s in samples}
        # patient 级 OS 临床
        clin = http_get_json(f"{CBIO}/studies/{study}/clinical-data?clinicalDataType=PATIENT&pageSize=10000")
        osd = {}
        for rec in clin:
            attr = rec.get("clinicalAttributeId")
            if attr in ("OS_MONTHS", "OS_STATUS"):
                osd.setdefault(rec["patientId"], {})[attr] = rec.get("value")
        # patient -> 样本(取 01)
        rows = {}
        for pid, d in osd.items():
            sid = next((s for s, p in sp.items() if p == pid), None)
            if sid and "OS_MONTHS" in d and "OS_STATUS" in d:
                rows[sid] = {"OS_MONTHS": d["OS_MONTHS"], "OS_STATUS": d["OS_STATUS"]}
        df = pd.DataFrame.from_dict(rows, orient="index")
        p = os.path.join(BASE, f"TCGA.{cancer}.sampleMap", f"{cancer}_clinicalMatrix.gz")
        df.to_csv(p, sep="\t", compression="gzip")
        print(f"  写出 {p} ({len(df)} 样本)")

if __name__ == "__main__":
    print("=== 1/2 构建表达矩阵（本地真实 TPM）===")
    build_expression()
    print("\n=== 2/2 构建临床矩阵（cBioPortal 真实 OS）===")
    build_clinical()
    print("\n适配输入就绪。可运行 h2_h3_tcga.py 端到端验证。")
