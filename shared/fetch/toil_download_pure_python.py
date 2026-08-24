#!/usr/bin/env python3
"""
Pure Python streaming download of ZP3 transcript isoform TPM (Toil Hub)

Avoid shell pipes, filter directly with requests + gzip stream
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
    "ENST00000394860": "ZP3_CANCER_CANDIDATE_5exon_8593bp",  # shortest protein-coding, candidate lacking signal peptide
    "ENST00000466960": "ZP3_retained_intron",
    "ENST00000479793": "ZP3_CDS_not_defined_2exon",
    "ENST00000467555": "ZP3_CDS_not_defined_4exon",
    "ENST00001135277": "ZP3_transcript_8exon_17148bp",
}

def main():
    print("=" * 60)
    print("Pure Python streaming download: ZP3 transcript isoform TPM")
    print(f"URL: {TOIL_URL}")
    print(f"Target transcripts: {len(ZP3_TRANSCRIPTS)}")
    print("=" * 60)
    
    # Streaming download + line-by-line filtering
    print("\n[1/3] Streaming download + line-by-line filtering of ZP3 transcripts...")
    t0 = time.time()
    
    header = None
    found_lines = {}
    total_bytes = 0
    last_report = 0
    
    enst_set = set(ZP3_TRANSCRIPTS.keys())
    
    with requests.get(TOIL_URL, stream=True, timeout=(30, 600)) as r:
        r.raise_for_status()
        content_length = int(r.headers.get('Content-Length', 0))
        print(f"  File size: {content_length / 1e9:.2f} GB")
        print("  Starting streaming read (estimated 10-40 minutes)...")
        
        # Use iterator to read raw bytes
        raw_iter = r.iter_content(chunk_size=65536)  # 64KB chunks
        
        # Use zlib.decompressobj() for streaming gzip decompression
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
            
            # Process complete lines in buffer
            while b"\n" in buffer:
                line_bytes, buffer = buffer.split(b"\n", 1)
                total_lines += 1
                
                try:
                    line = line_bytes.decode("utf-8", errors="replace")
                except:
                    continue
                
                if total_lines == 1:
                    # First line is the header
                    header = line
                    print(f"  Header: {line[:100]}...")
                    # Save header
                    with open(OUTPUT_HEADER, "w") as fh:
                        fh.write(line + "\n")
                    continue
                
                # Check if it is a ZP3 transcript
                for tid in enst_set:
                    if line.startswith(tid + "\t"):
                        found_lines[tid] = line
                        found_count += 1
                        break
            
            # Progress
            if total_bytes - last_report > 50_000_000:  # every 50MB
                elapsed = time.time() - t0
                pct = total_bytes / content_length * 100 if content_length else 0
                mb = total_bytes / 1e6
                rate = mb / (elapsed / 60) if elapsed > 0 else 0
                print(f"  [{elapsed/60:.1f} min] {mb:.0f} MB / {content_length/1e6:.0f} MB ({pct:.1f}%), ~{rate:.1f} MB/min, found {found_count} lines")
                last_report = total_bytes
        
        # Process the remaining buffer
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
    print(f"  Streaming read complete! Elapsed: {elapsed/60:.1f} min, total lines: {total_lines}, found: {found_count}")
    
    # Step 2: Write results
    print("\n[2/3] Writing result file...")
    
    with open(OUTPUT_TSV, "w", encoding="utf-8") as f:
        # Write header
        if header:
            f.write(header + "\n")
        # Write ZP3 rows
        for tid in sorted(found_lines.keys()):
            f.write(found_lines[tid] + "\n")
    
    size_on_disk = os.path.getsize(OUTPUT_TSV)
    print(f"  Output file: {OUTPUT_TSV}")
    print(f"  Size: {size_on_disk / 1e6:.2f} MB")
    
    # Report found transcripts
    print(f"\n  Found transcripts ({len(found_lines)}/{len(ZP3_TRANSCRIPTS)}):")
    for tid in sorted(found_lines.keys()):
        print(f"    ✓ {tid} ({ZP3_TRANSCRIPTS[tid]})")
    for tid in sorted(set(ZP3_TRANSCRIPTS) - set(found_lines.keys())):
        print(f"    ✗ {tid} ({ZP3_TRANSCRIPTS[tid]}) — not found in data (expression may be 0)")
    
    # Step 3: Metadata
    print("\n[3/3] Generating metadata...")
    with open(OUTPUT_META, "w", encoding="utf-8") as f:
        f.write(f"# ZP3 transcript isoform TPM data\n")
        f.write(f"# Source: {TOIL_URL}\n")
        f.write(f"# Download time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Elapsed: {elapsed/60:.1f} minutes\n")
        f.write(f"# Total lines: {total_lines}\n")
        f.write(f"# Transcripts found: {len(found_lines)}/{len(ZP3_TRANSCRIPTS)}\n")
        f.write(f"# Data format: log2(TPM + 0.001) RSEM, column 1=transcript_id, subsequent columns=samples\n")
        f.write(f"# Transcript list:\n")
        for tid in sorted(found_lines.keys()):
            f.write(f"#   {tid} = {ZP3_TRANSCRIPTS[tid]}\n")
        for tid in sorted(set(ZP3_TRANSCRIPTS) - set(found_lines.keys())):
            f.write(f"#   (missing) {tid} = {ZP3_TRANSCRIPTS[tid]}\n")
    
    print(f"  Metadata: {OUTPUT_META}")
    print("\n✓ Done!")
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
