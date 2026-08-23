#!/usr/bin/env python3
"""
UCSC Xena Toil Hub: 流式下载 + 过滤 ZP3 转录本 isoform TPM 数据

Toil Hub 数据集: TcgaTargetGtex_rsem_isoform_tpm (4.47 GB)
包含所有 TCGA + GTEx 样本的 RSEM 转录本 TPM (log2(TPM+0.001))

策略: 流式 gzip 读取, 只保留 ZP3 基因 (ENSG00000188372) 的转录本行

ZP3 转录本 (Ensembl GRCh38):
- ENST00000336517 (protein_coding, 9 exon) — ZP3-202?
- ENST00000394857 (protein_coding, 8 exon, CANONICAL) — 经典分泌型
- ENST00000416245 (protein_coding, 7 exon)
- ENST00000394860 (protein_coding, 5 exon) — 最短蛋白编码, 可能缺信号肽 = ZP3-Cancer 候选
- ENST00000466960 (retained_intron)
- ENST00000479793 (protein_coding_CDS_not_defined, 2 exon)
- ENST00000467555 (protein_coding_CDS_not_defined, 4 exon)
- ENST00001135277 (protein_coding, 8 exon)
"""

import gzip
import subprocess
import sys
import os
import time
import csv

TOIL_URL = "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_isoform_tpm.gz"
OUTPUT_TSV = os.path.join(os.path.dirname(__file__), "zp3_toil_isoform_tpm.tsv")
OUTPUT_META = os.path.join(os.path.dirname(__file__), "zp3_toil_meta.txt")
SPEED_SAMPLE_MB = 10  # 采样多少 MB 估算速度

# ZP3 转录本 ID
ZP3_TRANSCRIPTS = {
    "ENST00000336517": "ZP3_transcript_9exon_44574bp",
    "ENST00000394857": "ZP3_CANONICAL_8exon_17119bp",
    "ENST00000416245": "ZP3_transcript_7exon_13433bp",
    "ENST00000394860": "ZP3_CANCER_CANDIDATE_5exon_8593bp",
    "ENST00000466960": "ZP3_retained_intron",
    "ENST00000479793": "ZP3_CDS_not_defined_2exon",
    "ENST00000467555": "ZP3_CDS_not_defined_4exon",
    "ENST00001135277": "ZP3_transcript_8exon_17148bp",
}

def main():
    print("=" * 60)
    print("Toil Hub 流式下载: ZP3 转录本 isoform TPM")
    print(f"URL: {TOIL_URL}")
    print(f"ZP3 转录本数: {len(ZP3_TRANSCRIPTS)}")
    print("=" * 60)
    
    # 步骤1: 下载并过滤
    print("\n[1/3] 流式下载 + 过滤 ZP3 转录本...")
    t0 = time.time()
    
    # 构建 grep 正则: 匹配以任一个 ENST ID 开头的行
    enst_ids = "|".join(ZP3_TRANSCRIPTS.keys())
    grep_pattern = f"^{enst_ids}\\t"
    
    # 使用 curl 管道 zcat 管道 grep
    cmd = f'curl -s -L "{TOIL_URL}" | zcat 2>/dev/null | grep -E "{grep_pattern}" > "{OUTPUT_TSV}"'
    
    print(f"  命令: curl ... | zcat | grep -E '{grep_pattern}'")
    print("  注: 需流式读取完整 4.47 GB 文件...")
    
    process = subprocess.Popen(cmd, shell=True, executable="/bin/bash",
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # 等待完成 (不设超时)
    print("  等待下载完成 (可能需要 10-60 分钟)...")
    
    # 定期打印进度
    last_size = 0
    stall_count = 0
    while process.poll() is None:
        time.sleep(15)
        if os.path.exists(OUTPUT_TSV):
            current_size = os.path.getsize(OUTPUT_TSV)
            if current_size > last_size:
                elapsed = time.time() - t0
                mb = current_size / 1e6
                rate = mb / (elapsed / 60) if elapsed > 0 else 0
                print(f"  [{(elapsed/60):.1f} min] 已读取 {mb:.2f} MB, 速率 ~{rate:.1f} MB/min")
                last_size = current_size
                stall_count = 0
            else:
                stall_count += 1
        else:
            stall_count += 1
        
        if stall_count > 40:  # 10 分钟无进度
            print("  ⚠ 长时间无进度, 继续等待...")
            stall_count = 0
    
    elapsed = time.time() - t0
    print(f"  下载完成! 耗时: {elapsed/60:.1f} 分钟")
    
    # 步骤2: 检查结果
    print("\n[2/3] 检查结果...")
    if os.path.exists(OUTPUT_TSV):
        size_mb = os.path.getsize(OUTPUT_TSV) / 1e6
        print(f"  输出文件大小: {size_mb:.2f} MB")
        
        with open(OUTPUT_TSV) as f:
            lines = f.readlines()
        print(f"  总行数: {len(lines)}")
        
        found_transcripts = set()
        for line in lines:
            for tid in ZP3_TRANSCRIPTS:
                if line.startswith(tid):
                    found_transcripts.add(tid)
        
        print(f"  找到转录本: {len(found_transcripts)}/{len(ZP3_TRANSCRIPTS)}")
        for tid in sorted(found_transcripts):
            print(f"    ✓ {tid} ({ZP3_TRANSCRIPTS[tid]})")
        
        missing = set(ZP3_TRANSCRIPTS) - found_transcripts
        for tid in sorted(missing):
            print(f"    ✗ {tid} ({ZP3_TRANSCRIPTS[tid]}) — 未在数据中找到")
        
        # 步骤3: 生成元数据
        print("\n[3/3] 生成元数据...")
        with open(OUTPUT_META, "w") as f:
            f.write(f"# ZP3 转录本 isoform TPM 数据\n")
            f.write(f"# 来源: {TOIL_URL}\n")
            f.write(f"# 下载时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# 耗时: {elapsed/60:.1f} 分钟\n")
            f.write(f"# 行数: {len(lines)}\n")
            f.write(f"# 找到转录本: {len(found_transcripts)}/{len(ZP3_TRANSCRIPTS)}\n")
            f.write(f"# 转录本列表:\n")
            for tid in sorted(found_transcripts):
                f.write(f"#   {tid} = {ZP3_TRANSCRIPTS[tid]}\n")
        
        print(f"  完成! 产物: {OUTPUT_TSV}, {OUTPUT_META}")
        print(f"  文件大小: {size_mb:.2f} MB")
        print(f"  行数: {len(lines)}")
        
        return True
    else:
        print("  ⚠ 输出文件未生成, 下载可能失败")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
