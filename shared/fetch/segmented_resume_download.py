#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分段续传下载器 (Segmented Resume Downloader)
============================================
解决后台进程被沙箱回收导致大文件下载中断的问题：
- 每个下载任务按 SEG 大小分 N 段，每段独立下载到 .part{N} 文件
- 进程被回收后重启本脚本，自动检测已完成的 .part 文件并跳过，只补未完成的段
- 所有段完成后合并为最终文件，校验大小，写 *_DONE.txt 标记

支持 HTTP Range (Accept-Ranges: bytes)。已验证 toil.xenahubs.net 与
datasets.cellxgene.cziscience.com 均支持。

用法: python segmented_resume_download.py
"""
import os
import sys
import time
import requests

BASE = os.path.dirname(os.path.abspath(__file__))
SEG = 400 * 1024 * 1024  # 每段 400 MB
MAX_RETRY = 3            # 每段最大重试次数

# 任务定义： (名称, URL, 输出相对路径, 断点续传文件(可选))
# 顺序按优先级：GTEx(文章1必需) -> 黑色素瘤(跨癌种验证) -> Toil(文章3)
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
    """带重试的 Range 下载单段"""
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
            log(f'    段下载异常(尝试{attempt}/{MAX_RETRY}): {str(e)[:80]}')
            if os.path.exists(part):
                os.remove(part)
            time.sleep(5)
    raise last_err


def process_task(name, url, output_rel, resume_rel):
    output = os.path.join(BASE, output_rel)
    os.makedirs(os.path.dirname(output), exist_ok=True)

    # 断点续传文件处理：把已有文件重命名为 .resume（避免与 output 路径冲突）
    resume_file = None
    if resume_rel:
        resume_path = os.path.join(BASE, resume_rel)
        if os.path.exists(resume_path) and not os.path.exists(output):
            resume_file = resume_path
        elif os.path.exists(output) and not os.path.exists(resume_path):
            # output 已存在但 resume 不存在：说明之前合并未完成或异常，重命名为 resume
            os.rename(output, resume_path)
            resume_file = resume_path

    total = get_size(url)
    log(f'[{name}] 总大小 {total/1e9:.2f} GB')

    start_offset = 0
    if resume_file and os.path.exists(resume_file):
        start_offset = os.path.getsize(resume_file)
        log(f'[{name}] 检测到断点文件 {start_offset/1e6:.0f} MB，从该位置续传')

    # 计算分段
    segments = []
    pos = start_offset
    idx = 0
    while pos < total:
        end = min(pos + SEG - 1, total - 1)
        segments.append((idx, pos, end))
        pos = end + 1
        idx += 1
    log(f'[{name}] 分 {len(segments)} 段 (每段 {SEG/1e6:.0f} MB)')

    for seg_idx, s, e in segments:
        part = f'{output}.part{seg_idx}'
        expected = e - s + 1
        if os.path.exists(part) and os.path.getsize(part) == expected:
            log(f'  [{name}] 段{seg_idx} 已完成({expected/1e6:.0f}MB)，跳过')
            continue
        if os.path.exists(part):
            os.remove(part)
        log(f'  [{name}] 段{seg_idx}: bytes {s}-{e} ({expected/1e6:.0f}MB) 下载中...')
        t0 = time.time()
        download_range(url, s, e, part)
        dt = time.time() - t0
        speed = (expected / 1e6) / dt if dt > 0 else 0
        log(f'  [{name}] 段{seg_idx} 完成 ({speed:.0f} MB/s)')

    # 合并
    log(f'[{name}] 合并 {len(segments)} 段...')
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
        log(f'[{name}] ✅ 完成，大小校验一致 {actual} bytes')
        with open(os.path.join(BASE, f'{name}_DONE.txt'), 'w') as f:
            f.write(str(actual))
    else:
        log(f'[{name}] ❌ 大小不匹配! 实际 {actual} vs 期望 {total}')


def main():
    log('=== 分段续传下载器启动 (pid=%d) ===' % os.getpid())
    for name, url, out_rel, resume_rel in TASKS:
        done_mark = os.path.join(BASE, f'{name}_DONE.txt')
        if os.path.exists(done_mark):
            log(f'[{name}] 已有完成标记，跳过')
            continue
        # 任务级重试：应对瞬时网络错误，进程被杀则靠外部重启续传
        for attempt in range(1, 4):
            try:
                process_task(name, url, out_rel, resume_rel)
                break
            except Exception as e:
                log(f'[{name}] ⚠️ 任务异常(尝试{attempt}/3): {str(e)[:120]}')
                if attempt < 3:
                    log(f'[{name}] 30秒后重试...')
                    time.sleep(30)
                else:
                    log(f'[{name}] 本轮放弃，下次重启将自动续传未完成段')
    log('=== 所有任务处理完毕 ===')
    with open(os.path.join(BASE, 'ALL_DOWNLOADED.txt'), 'w') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M:%S'))


if __name__ == '__main__':
    main()
