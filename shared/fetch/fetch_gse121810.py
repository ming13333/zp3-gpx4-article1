# -*- coding: utf-8 -*-
import os as _os
def _project_root():
    d = _os.path.dirname(_os.path.abspath(__file__))
    while True:
        if _os.path.isdir(_os.path.join(d, "output")):
            return d
        p = _os.path.dirname(d)
        if p == d:
            break
        d = p
    return _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
ROOT = _project_root()
"""
获取 GSE121810 胶质瘤免疫治疗数据集
数据来源：GEO FTP 补充文件
"""
import os
import requests
import pandas as pd
import numpy as np

OUT_DIR = os.path.join(ROOT, "output", "h2_bulk")
os.makedirs(OUT_DIR, exist_ok=True)

# 探测补充文件
print('探测 GSE121810 补充文件...')
url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121810/suppl/'
try:
    r = requests.get(url, timeout=30)
    print(f'状态码: {r.status_code}')
    # 解析 HTML 列表
    import re
    files = re.findall(r'href="([^"]+)"', r.text)
    print(f'找到 {len(files)} 个文件:')
    for f in files:
        print(f'  {f}')
except Exception as e:
    print(f'探测失败: {e}')
    exit(1)

# 下载补充文件
print('\n下载补充文件...')
for file in files:
    if file.endswith('.xlsx') or file.endswith('.csv') or file.endswith('.txt.gz'):
        file_url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121810/suppl/{file}'
        local_file = os.path.join(OUT_DIR, file)
        print(f'下载 {file}...')
        try:
            r = requests.get(file_url, timeout=300)
            r.raise_for_status()
            with open(local_file, 'wb') as f:
                f.write(r.content)
            print(f'  保存: {local_file} ({len(r.content)} bytes)')
        except Exception as e:
            print(f'  下载失败: {e}')

# 读取下载的文件
print('\n读取下载的文件...')
for file in os.listdir(OUT_DIR):
    if file.startswith('GSE121810'):
        file_path = os.path.join(OUT_DIR, file)
        print(f'\n=== {file} ===')
        try:
            if file.endswith('.xlsx'):
                df = pd.read_excel(file_path, nrows=5)
                print(f'形状: {df.shape}')
                print(f'列名: {list(df.columns)[:10]}...')
                print(df.head(2).to_string())
            elif file.endswith('.csv'):
                df = pd.read_csv(file_path, nrows=5)
                print(f'形状: {df.shape}')
                print(f'列名: {list(df.columns)[:10]}...')
                print(df.head(2).to_string())
            elif file.endswith('.txt.gz'):
                df = pd.read_csv(file_path, sep='\t', nrows=5, compression='gzip')
                print(f'形状: {df.shape}')
                print(f'列名: {list(df.columns)[:10]}...')
                print(df.head(2).to_string())
        except Exception as e:
            print(f'读取失败: {e}')

# 下载 series matrix
print('\n下载 series matrix...')
series_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121810/matrix/GSE121810_series_matrix.txt.gz'
series_file = os.path.join(OUT_DIR, 'GSE121810_series_matrix.txt.gz')
try:
    r = requests.get(series_url, timeout=300)
    r.raise_for_status()
    with open(series_file, 'wb') as f:
        f.write(r.content)
    print(f'保存: {series_file} ({len(r.content)} bytes)')
except Exception as e:
    print(f'下载失败: {e}')