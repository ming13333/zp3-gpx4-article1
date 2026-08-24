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
Fetch GSE121810 glioma immunotherapy dataset
Data source: GEO FTP supplementary files
"""
import os
import requests
import pandas as pd
import numpy as np

OUT_DIR = os.path.join(ROOT, "output", "h2_bulk")
os.makedirs(OUT_DIR, exist_ok=True)

# Probe supplementary files
print('Probing GSE121810 supplementary files...')
url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121810/suppl/'
try:
    r = requests.get(url, timeout=30)
    print(f'Status code: {r.status_code}')
    # Parse HTML list
    import re
    files = re.findall(r'href="([^"]+)"', r.text)
    print(f'Found {len(files)} files:')
    for f in files:
        print(f'  {f}')
except Exception as e:
    print(f'Probe failed: {e}')
    exit(1)

# Download supplementary files
print('\nDownloading supplementary files...')
for file in files:
    if file.endswith('.xlsx') or file.endswith('.csv') or file.endswith('.txt.gz'):
        file_url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121810/suppl/{file}'
        local_file = os.path.join(OUT_DIR, file)
        print(f'Downloading {file}...')
        try:
            r = requests.get(file_url, timeout=300)
            r.raise_for_status()
            with open(local_file, 'wb') as f:
                f.write(r.content)
            print(f'  Saved: {local_file} ({len(r.content)} bytes)')
        except Exception as e:
            print(f'  Download failed: {e}')

# Read downloaded files
print('\nReading downloaded files...')
for file in os.listdir(OUT_DIR):
    if file.startswith('GSE121810'):
        file_path = os.path.join(OUT_DIR, file)
        print(f'\n=== {file} ===')
        try:
            if file.endswith('.xlsx'):
                df = pd.read_excel(file_path, nrows=5)
                print(f'Shape: {df.shape}')
                print(f'Columns: {list(df.columns)[:10]}...')
                print(df.head(2).to_string())
            elif file.endswith('.csv'):
                df = pd.read_csv(file_path, nrows=5)
                print(f'Shape: {df.shape}')
                print(f'Columns: {list(df.columns)[:10]}...')
                print(df.head(2).to_string())
            elif file.endswith('.txt.gz'):
                df = pd.read_csv(file_path, sep='\t', nrows=5, compression='gzip')
                print(f'Shape: {df.shape}')
                print(f'Columns: {list(df.columns)[:10]}...')
                print(df.head(2).to_string())
        except Exception as e:
            print(f'Read failed: {e}')

# Download series matrix
print('\nDownloading series matrix...')
series_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121810/matrix/GSE121810_series_matrix.txt.gz'
series_file = os.path.join(OUT_DIR, 'GSE121810_series_matrix.txt.gz')
try:
    r = requests.get(series_url, timeout=300)
    r.raise_for_status()
    with open(series_file, 'wb') as f:
        f.write(r.content)
    print(f'Saved: {series_file} ({len(r.content)} bytes)')
except Exception as e:
    print(f'Download failed: {e}')
