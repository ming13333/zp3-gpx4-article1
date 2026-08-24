#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Segmented Resume Downloader (Segmented Resume Downloader)
============================================
Solve the problem of large file downloads being interrupted when background processes are reclaimed by the sandbox:
- Each download task is divided into N segments by SEG size, and each segment is downloaded independently to .part{N} file
- After the process is reclaimed, restart this script to automatically detect completed .part files and skip them, only downloading missing segments
- After all segments are completed, merge into the final file, verify the size, and write a *_DONE.txt marker

Supports HTTP Range (Accept-Ranges: bytes). Verified toil.xenahubs.net and
datasets.cellxgene.cziscience.com both support it.

Usage: python segmented_resume_download.py
"""
import os
import sys
import time
import requests

BASE = os.path.dirname(os.path.abspath(__file__))
# Each segment is 400 MB
# Maximum retry count per segment

# Task definition: (name, URL, output relative path, resume file (optional))
# Order by priority: GTEx(Article 1 required) -> melanoma(cross-cancer validation) -> Toil(Article 3)
TASKS = [
    ("GTEx_gene",
     "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz",
     "TcgaTargetGtex_rsem_gene_tpm.gz", None),
    ("GTEx_pheno",
     "https://toil.xenahubs.net/download/GTEX_phenotype.gz",
     "GTEX_phenotype.gz", None),
    ("Melanoma",
     "https://datasets.cellxgene.cziscience.com/1b76227b-c731-4807-9487-ad5e4d24e0d0.h5ad",
     "sc_data/melanoma_myeloid.h5ad", "sc_data/melanoma_myeloid.h5ad.resume"),
    ("Toil_isoform",
     "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_isoform_tpm.gz",
     "TcgaTargetGtex_rsem_isoform_tpm.gz", None),
]


def log(msg):
    ts = time.strftime('%Y-%m-%d %H:%M:%S')
    line = f'[{ts}] {msg}'
    print(line, flush=True)
    with open(os.path.join(BASE, 'download_progress.log'), 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def get_size(url):
    r = requests.head(url, timeout=30, allow_redirects=True)
    r.raise_for_status()
    return int(r.headers['Content-Length'])


def download_range(url, start, end, part):
    """Download a single range with retry"""
    headers = {'Range': f'bytes={start}-{end}'}
    last_err = None
    for attempt in range(1, MAX_RETRY + 1):
        try:
            r = requests.get(url, headers=headers, stream=True, timeout=(30, 600))
            r.raise_for_status()
            with open(part, 'wb') as f:
                for chunk in r.iter_content(1024 * 1024):
                    if chunk:
                        f.write(chunk)
            return True
        except Exception as e:
            last_err = e
            log(f'    Segment download error (attempt {attempt}/{MAX_RETRY}): {str(e)[:80]}')
            if os.path.exists(part):
                os.remove(part)
            time.sleep(5)
    raise last_err


def process_task(name, url, output_rel, resume_rel):
    output = os.path.join(BASE, output_rel)
    os.makedirs(os.path.dirname(output), exist_ok=True)

    # Resume file handling: rename existing file to .resume (avoid conflict with output path)
    resume_file = None
    if resume_rel:
        resume_path = os.path.join(BASE, resume_rel)
        if os.path.exists(resume_path) and not os.path.exists(output):
            resume_file = resume_path
        elif os.path.exists(output) and not os.path.exists(resume_path):
            # output exists but resume doesn't: means previous merge was incomplete or abnormal, rename to resume
            os.rename(output, resume_path)
            resume_file = resume_path

    total = get_size(url)
    log(f'[{name}] Total size {total/1e9:.2f} GB')

    start_offset = 0
    if resume_file and os.path.exists(resume_file):
        start_offset = os.path.getsize(resume_file)
        log(f'[{name}] Detected resume file {start_offset/1e6:.0f} MB, resuming from that position')

    # Calculate segments
    segments = []
    pos = start_offset
    idx = 0
    while pos < total:
        end = min(pos + SEG - 1, total - 1)
        segments.append((idx, pos, end))
        pos = end + 1
        idx += 1
    log(f'[{name}] Split into {len(segments)} segments (each {SEG/1e6:.0f} MB)')

    for seg_idx, s, e in segments:
        part = f'{output}.part{seg_idx}'
        expected = e - s + 1
        if os.path.exists(part) and os.path.getsize(part) == expected:
            log(f'  [{name}] Segment {seg_idx} already done ({expected/1e6:.0f}MB), skipping')
            continue
        if os.path.exists(part):
            os.remove(part)
        log(f'  [{name}] Segment {seg_idx}: bytes {s}-{e} ({expected/1e6:.0f}MB) downloading...')
        t0 = time.time()
        download_range(url, s, e, part)
        dt = time.time() - t0
        speed = (expected / 1e6) / dt if dt > 0 else 0
        log(f'  [{name}] Segment {seg_idx} done ({speed:.0f} MB/s)')

    # Merge
    log(f'[{name}] Merging {len(segments)} segments...')
    with open(output, 'wb') as out:
        if resume_file and os.path.exists(resume_file):
            with open(resume_file, 'rb') as rf:
                while True:
                    chunk = rf.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            os.remove(resume_file)
        for seg_idx, s, e in segments:
            part = f'{output}.part{seg_idx}'
            with open(part, 'rb') as pf:
                while True:
                    chunk = pf.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            os.remove(part)

    actual = os.path.getsize(output)
    if actual == total:
        log(f'[{name}] ✅ Done, size verification consistent {actual} bytes')
        with open(os.path.join(BASE, f'{name}_DONE.txt'), 'w') as f:
            f.write(str(actual))
    else:
        log(f'[{name}] ❌ Size mismatch! Actual {actual} vs expected {total}')


def main():
    log('=== Segmented resume downloader started (pid=%d) ===' % os.getpid())
    for name, url, out_rel, resume_rel in TASKS:
        done_mark = os.path.join(BASE, f'{name}_DONE.txt')
        if os.path.exists(done_mark):
            log(f'[{name}] Already has completion marker, skipping')
            continue
        # Task-level retry: handle transient network errors; if the process is killed, external restart resumes the transfer
        for attempt in range(1, 4):
            try:
                process_task(name, url, out_rel, resume_rel)
                break
            except Exception as e:
                log(f'[{name}] ⚠️ Task exception (attempt {attempt}/3): {str(e)[:120]}')
                if attempt < 3:
                    log(f'[{name}] Retrying in 30 seconds...')
                    time.sleep(30)
                else:
                    log(f'[{name}] Giving up this round, next restart will automatically resume incomplete segments')
    log('=== All tasks processed ===')
    with open(os.path.join(BASE, 'ALL_DOWNLOADED.txt'), 'w') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M:%S'))


if __name__ == '__main__':
    main()
