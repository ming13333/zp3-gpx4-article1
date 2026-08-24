#!/usr/bin/env python3
"""
UCSC Xena Toil Hub: streaming download + filter ZP3 transcript isoform TPM data

Toil Hub dataset: TcgaTargetGtex_rsem_isoform_tpm (4.47 GB)
Contains RSEM transcript TPM for all TCGA + GTEx samples (log2(TPM+0.001))

Strategy: stream gzip read, keep only transcript rows for ZP3 gene (ENSG00000188372)

ZP3 transcripts (Ensembl GRCh38):
- ENST00000336517 (protein_coding, 9 exon) — ZP3-202?
- ENST00000394857 (protein_coding, 8 exon, CANONICAL) — canonical secreted form
- ENST00000416245 (protein_coding, 7 exon)
- ENST00000394860 (protein_coding, 5 exon) — shortest protein-coding, may lack signal peptide = ZP3-Cancer candidate
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
SPEED_SAMPLE_MB = 10  # sample how many MB to estimate speed

# ZP3 transcript IDs
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
    print("Toil Hub streaming download: ZP3 transcript isoform TPM")
    print(f"URL: {TOIL_URL}")
    print(f"ZP3 transcript count: {len(ZP3_TRANSCRIPTS)}")
    print("=" * 60)
    
    # Step 1: download and filter
    print("\n[1/3] Streaming download + filter ZP3 transcripts...")
    t0 = time.time()
    
    # Build grep regex: match lines starting with any ENST ID
    enst_ids = "|".join(ZP3_TRANSCRIPTS.keys())
    grep_pattern = f"^{enst_ids}\\t"
    
    # Use curl piped to zcat piped to grep
    cmd = f'curl -s -L "{TOIL_URL}" | zcat 2>/dev/null | grep -E "{grep_pattern}" > "{OUTPUT_TSV}"'
    
    print(f"  Command: curl ... | zcat | grep -E '{grep_pattern}'")
    print("  Note: needs to stream-read the full 4.47 GB file...")
    
    process = subprocess.Popen(cmd, shell=True, executable="/bin/bash",
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Wait for completion (no timeout)
    print("  Waiting for download to complete (may take 10-60 minutes)...")
    
    # Periodically print progress
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
                print(f"  [{(elapsed/60):.1f} min] Read {mb:.2f} MB, rate ~{rate:.1f} MB/min")
                last_size = current_size
                stall_count = 0
            else:
                stall_count += 1
        else:
            stall_count += 1
        
        if stall_count > 40:  # 10 minutes without progress
            print("  ⚠ No progress for a long time, keep waiting...")
            stall_count = 0
    
    elapsed = time.time() - t0
    print(f"  Download complete! Time spent: {elapsed/60:.1f} minutes")
    
    # Step 2: check results
    print("\n[2/3] Checking results...")
    if os.path.exists(OUTPUT_TSV):
        size_mb = os.path.getsize(OUTPUT_TSV) / 1e6
        print(f"  Output file size: {size_mb:.2f} MB")
        
        with open(OUTPUT_TSV) as f:
            lines = f.readlines()
        print(f"  Total lines: {len(lines)}")
        
        found_transcripts = set()
        for line in lines:
            for tid in ZP3_TRANSCRIPTS:
                if line.startswith(tid):
                    found_transcripts.add(tid)
        
        print(f"  Found transcripts: {len(found_transcripts)}/{len(ZP3_TRANSCRIPTS)}")
        for tid in sorted(found_transcripts):
            print(f"    ✓ {tid} ({ZP3_TRANSCRIPTS[tid]})")
        
        missing = set(ZP3_TRANSCRIPTS) - found_transcripts
        for tid in sorted(missing):
            print(f"    ✗ {tid} ({ZP3_TRANSCRIPTS[tid]}) — not found in data")
        
        # Step 3: Generate metadata
        print("\n[3/3] Generating metadata...")
        with open(OUTPUT_META, "w") as f:
            f.write(f"# ZP3 transcript isoform TPM data\n")
            f.write(f"# Source: {TOIL_URL}\n")
            f.write(f"# Download time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Time elapsed: {elapsed/60:.1f} minutes\n")
            f.write(f"# Line count: {len(lines)}\n")
            f.write(f"# Transcripts found: {len(found_transcripts)}/{len(ZP3_TRANSCRIPTS)}\n")
            f.write(f"# Transcript list:\n")
            for tid in sorted(found_transcripts):
                f.write(f"#   {tid} = {ZP3_TRANSCRIPTS[tid]}\n")
        
        print(f"  Done! Output files: {OUTPUT_TSV}, {OUTPUT_META}")
        print(f"  File size: {size_mb:.2f} MB")
        print(f"  Lines: {len(lines)}")
        
        return True
    else:
        print("  ⚠ Output file not generated, download may have failed")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
