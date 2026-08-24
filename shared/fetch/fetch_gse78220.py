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
Fetch GSE78220 glioma immunotherapy dataset
Data source: GEO FTP supplementary file GSE78220_PatientFPKM.xlsx
"""
import os
import requests
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# Set paths
OUT_DIR = os.path.join(ROOT, "output", "h2_bulk")
os.makedirs(OUT_DIR, exist_ok=True)

# Download file
url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx'
local_file = os.path.join(OUT_DIR, 'GSE78220_PatientFPKM.xlsx')

print('Downloading GSE78220 data...')
try:
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    with open(local_file, 'wb') as f:
        f.write(r.content)
    print(f'Download complete: {local_file} ({len(r.content)} bytes)')
except Exception as e:
    print(f'Download failed: {e}')
    exit(1)

# Read Excel file
print('\nReading Excel file...')
try:
    # Try to read all sheets
    xl = pd.ExcelFile(local_file)
    print(f'Sheet names: {xl.sheet_names}')
    
    # Usually the first sheet is the expression matrix, the second is clinical data
    # First look at the first few rows of each sheet
    for sheet in xl.sheet_names:
        df = pd.read_excel(local_file, sheet_name=sheet, nrows=5)
        print(f'\n=== Sheet: {sheet} ===')
        print(f'Shape: {df.shape}')
        print(f'Columns: {list(df.columns)[:10]}...')
        print(df.head(2).to_string())
        
except Exception as e:
    print(f'Read failed: {e}')
    exit(1)

# Assume the first sheet is the expression matrix, the second is clinical data
# Adjust according to the actual structure
print('\nParsing data...')

# Read the full data
try:
    # Expression matrix
    expr_df = pd.read_excel(local_file, sheet_name=0, index_col=0)
    print(f'Expression matrix shape: {expr_df.shape}')
    
    # Clinical data
    clin_df = pd.read_excel(local_file, sheet_name=1, index_col=0)
    print(f'Clinical data shape: {clin_df.shape}')
    print(f'Clinical data columns: {list(clin_df.columns)}')
    
    # Save as CSV
    expr_csv = os.path.join(OUT_DIR, 'gse78220_expression.csv')
    clin_csv = os.path.join(OUT_DIR, 'gse78220_clinical.csv')
    
    expr_df.to_csv(expr_csv)
    clin_df.to_csv(clin_csv)
    
    print(f'\nSaved:')
    print(f'Expression matrix: {expr_csv}')
    print(f'Clinical data: {clin_csv}')
    
    # Check ZP3 expression
    if 'ZP3' in expr_df.index:
        zp3_expr = expr_df.loc['ZP3']
        print(f'\nZP3 expression statistics:')
        print(f'Number of samples: {len(zp3_expr)}')
        print(f'Mean: {zp3_expr.mean():.3f}')
        print(f'Median: {zp3_expr.median():.3f}')
        print(f'Std: {zp3_expr.std():.3f}')
        print(f'Range: {zp3_expr.min():.3f} - {zp3_expr.max():.3f}')
    else:
        print('\nWarning: ZP3 gene not found')
        # Search for possible ZP3-related genes
        zp3_candidates = [g for g in expr_df.index if 'ZP3' in str(g).upper()]
        if zp3_candidates:
            print(f'Possible ZP3-related genes found: {zp3_candidates}')
    
    # Check treatment response information in clinical data
    print(f'\nClinical data preview:')
    print(clin_df.head(10).to_string())
    
    # Search for possible treatment response columns
    response_cols = [c for c in clin_df.columns if any(keyword in str(c).lower() for keyword in ['response', 'outcome', 'treatment', 'anti', 'pd', 'pd-l1', 'immunotherapy'])]
    if response_cols:
        print(f'\nPossible treatment response-related columns: {response_cols}')
        for col in response_cols[:3]:  # show only the first 3
            print(f'{col}: {clin_df[col].value_counts().to_dict()}')
    
    # Search for survival time columns
    survival_cols = [c for c in clin_df.columns if any(keyword in str(c).lower() for keyword in ['survival', 'os', 'overall', 'time', 'month', 'day', 'year'])]
    if survival_cols:
        print(f'\nPossible survival-related columns: {survival_cols}')
    
    # Search for status columns (death/survival)
    status_cols = [c for c in clin_df.columns if any(keyword in str(c).lower() for keyword in ['status', 'event', 'alive', 'dead', 'death', 'vital'])]
    if status_cols:
        print(f'\nPossible status-related columns: {status_cols}')
        
except Exception as e:
    print(f'Parsing failed: {e}')
    import traceback
    traceback.print_exc()
    exit(1)

print('\nData fetch complete.')
