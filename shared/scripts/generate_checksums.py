#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
生成输入数据校验和清单（I3 复现性修复的一部分）
================================================
为所有大型输入数据文件计算 sha256 与字节数，写入 CHECKSUMS.sha256。
二次运行的分析应当先校验本地文件与清单一致，方能保证可复现。
用法：
    python generate_checksums.py            # 依据下方 CURATED 列表生成
    python generate_checksums.py --verify   # 校验现有文件与清单是否一致
"""
import os
import sys
import hashlib
import json

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "shared", "data")

# 关键输入数据文件（相对 BASE）。这些文件由下载脚本真实获取，是分析的底层输入。
CURATED = [
    "phase1_knowledge_gap_filling/TcgaTargetGtex_rsem_gene_tpm.gz",
    "phase1_knowledge_gap_filling/TcgaTargetGtex_rsem_isoform_tpm.gz",
    "phase1_knowledge_gap_filling/GTEX_phenotype.gz",
    "cgga_validation/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",
    "cgga_validation/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",
    "cgga_validation/CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip",
    "cgga_validation/CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip",
    "h1_pilot/GSE84465_GBM_All_data.csv.gz",
    "h1_pilot/h1_adata.h5ad",
    "h1_pilot/h1_adata_subtyped.h5ad",
    "gse91061_validation/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",
    "immunotherapy_validation/exmat_censored_IMvigor210.csv",
    "phase1_knowledge_gap_filling/sc_data/rcc_htan.h5ad",
    "phase1_knowledge_gap_filling/sc_data/luad_htan.h5ad",
    "phase1_knowledge_gap_filling/sc_data/melanoma_myeloid.h5ad",
]


def sha256_of(path, chunk=8 * 1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main():
    mode = "--verify" if "--verify" in sys.argv else "--write"
    out = os.path.join(BASE, "CHECKSUMS.sha256")
    manifest = {}
    if mode == "--verify" and os.path.exists(out):
        with open(out) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                sha, rel = line.split(None, 1)
                manifest[rel] = sha

    records = []
    for rel in CURATED:
        fp = os.path.join(BASE, rel)
        if not os.path.exists(fp):
            print(f"  (跳过缺失) {rel}")
            continue
        print(f"  计算校验和: {rel} ...")
        sha = sha256_of(fp)
        size = os.path.getsize(fp)
        records.append((rel, sha, size))
        if mode == "--verify" and rel in manifest:
            ok = manifest[rel] == sha
            print(f"    [{'OK' if ok else 'MISMATCH'}] {rel}")

    if mode == "--write":
        with open(out, "w") as f:
            f.write("# sha256  <relative-path-from-output/>\n")
            f.write(f"# 生成环境: {sys.version.split()[0]}\n")
            for rel, sha, size in records:
                f.write(f"{sha}  {rel}  # {size} bytes\n")
        print(f"\n已写入 {len(records)} 条校验和 -> {out}")
    else:
        print(f"\n校验完成，共 {len(records)} 个文件。")


if __name__ == "__main__":
    main()
