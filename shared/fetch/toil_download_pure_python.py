#!/usr/bin/env python3
"""
纯 Python 流式下载 ZP3 转录本 isoform TPM (Toil Hub)

避免 shell pipe，直接用 requests + gzip stream 过滤
"""

import requests
import gzip
import io
import os
import time
import sys

TOIL_URL = "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_isoform_tpm.gz"
DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_TSV = os.path.join(DIR, "zp3_toil_isoform_tpm.tsv")
OUTPUT_META = os.path.join(DIR, "zp3_toil_meta.txt")
OUTPUT_HEADER = os.path.join(DIR, "zp3_toil_header.tsv")

ZP3_TRANSCRIPTS = {
    "ENST00000336517": "ZP3_transcript_9exon_44574bp",
    "ENST00000394857": "ZP3_CANONICAL_8exon_17119bp",
    "ENST00000416245": "ZP3_transcript_7exon_13433bp",
    "ENST00000394860": "ZP3_CANCER_CANDIDATE_5exon_8593bp",  # 最短蛋白编码, 缺信号肽候选
    "ENST00000466960": "ZP3_retained_intron",
    "ENST00000479793": "ZP3_CDS_not_defined_2exon",
    "ENST00000467555": "ZP3_CDS_not_defined_4exon",
    "ENST00001135277": "ZP3_transcript_8exon_17148bp",
}

def main():
    print("=" * 60)
    print("纯 Python 流式下载: ZP3 转录本 isoform TPM")
    print(f"URL: {TOIL_URL}")
    print(f"目标转录本: {len(ZP3_TRANSCRIPTS)}")
    print("=" * 60)
    
    # 流式下载 + 逐行过滤
    print("\n[1/3] 流式下载 + 逐行过滤 ZP3 转录本...")
    t0 = time.time()
    
    header = None
    found_lines = {}
    total_bytes = 0
    last_report = 0
    
    enst_set = set(ZP3_TRANSCRIPTS.keys())
    
    with requests.get(TOIL_URL, stream=True, timeout=(30, 600)) as r:
        r.raise_for_status()
        content_length = int(r.headers.get('Content-Length', 0))
        print(f"  文件大小: {content_length / 1e9:.2f} GB")
        print("  开始流式读取 (预计 10-40 分钟)...")
        
        # 使用迭代器读 raw bytes
        raw_iter = r.iter_content(chunk_size=65536)  # 64KB chunks
        
        # 使用 zlib.decompressobj() 做流式 gzip 解压
        import zlib
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)  # gzip mode
        
        buffer = b""
        line_buffer = b""
        total_lines = 0
        found_count = 0
        
        for chunk in raw_iter:
            total_bytes += len(chunk)
            decompressed = d.decompress(chunk)
            buffer += decompressed
            
            # 处理 buffer 中的完整行
            while b"\n" in buffer:
                line_bytes, buffer = buffer.split(b"\n", 1)
                total_lines += 1
                
                try:
                    line = line_bytes.decode("utf-8", errors="replace")
                except:
                    continue
                
                if total_lines == 1:
                    # 第一行是 header
                    header = line
                    print(f"  列头: {line[:100]}...")
                    # 保存 header
                    with open(OUTPUT_HEADER, "w") as fh:
                        fh.write(line + "\n")
                    continue
                
                # 检查是否为 ZP3 转录本
                for tid in enst_set:
                    if line.startswith(tid + "\t"):
                        found_lines[tid] = line
                        found_count += 1
                        break
            
            # 进度
            if total_bytes - last_report > 50_000_000:  # 每 50MB
                elapsed = time.time() - t0
                pct = total_bytes / content_length * 100 if content_length else 0
                mb = total_bytes / 1e6
                rate = mb / (elapsed / 60) if elapsed > 0 else 0
                print(f"  [{elapsed/60:.1f} min] {mb:.0f} MB / {content_length/1e6:.0f} MB ({pct:.1f}%), ~{rate:.1f} MB/min, 找到 {found_count} 行")
                last_report = total_bytes
        
        # 处理剩余的 buffer
        if buffer:
            total_lines += 1
            try:
                line = buffer.decode("utf-8", errors="replace")
                for tid in enst_set:
                    if line.startswith(tid + "\t"):
                        found_lines[tid] = line
                        found_count += 1
                        break
            except:
                pass
    
    elapsed = time.time() - t0
    print(f"  流式读取完成! 耗时: {elapsed/60:.1f} 分钟, 总行数: {total_lines}, 找到: {found_count}")
    
    # 步骤2: 写结果
    print("\n[2/3] 写入结果文件...")
    
    with open(OUTPUT_TSV, "w", encoding="utf-8") as f:
        # 写 header
        if header:
            f.write(header + "\n")
        # 写 ZP3 行
        for tid in sorted(found_lines.keys()):
            f.write(found_lines[tid] + "\n")
    
    size_on_disk = os.path.getsize(OUTPUT_TSV)
    print(f"  输出文件: {OUTPUT_TSV}")
    print(f"  大小: {size_on_disk / 1e6:.2f} MB")
    
    # 报告找到的转录本
    print(f"\n  找到转录本 ({len(found_lines)}/{len(ZP3_TRANSCRIPTS)}):")
    for tid in sorted(found_lines.keys()):
        print(f"    ✓ {tid} ({ZP3_TRANSCRIPTS[tid]})")
    for tid in sorted(set(ZP3_TRANSCRIPTS) - set(found_lines.keys())):
        print(f"    ✗ {tid} ({ZP3_TRANSCRIPTS[tid]}) — 未在数据中找到 (可能表达量为 0)")
    
    # 步骤3: 元数据
    print("\n[3/3] 生成元数据...")
    with open(OUTPUT_META, "w", encoding="utf-8") as f:
        f.write(f"# ZP3 转录本 isoform TPM 数据\n")
        f.write(f"# 来源: {TOIL_URL}\n")
        f.write(f"# 下载时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# 耗时: {elapsed/60:.1f} 分钟\n")
        f.write(f"# 总行数: {total_lines}\n")
        f.write(f"# 找到转录本: {len(found_lines)}/{len(ZP3_TRANSCRIPTS)}\n")
        f.write(f"# 数据格式: log2(TPM + 0.001) RSEM, 第1列=transcript_id, 后续列=样本\n")
        f.write(f"# 转录本列表:\n")
        for tid in sorted(found_lines.keys()):
            f.write(f"#   {tid} = {ZP3_TRANSCRIPTS[tid]}\n")
        for tid in sorted(set(ZP3_TRANSCRIPTS) - set(found_lines.keys())):
            f.write(f"#   (缺失) {tid} = {ZP3_TRANSCRIPTS[tid]}\n")
    
    print(f"  元数据: {OUTPUT_META}")
    print("\n✓ 完成!")
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
