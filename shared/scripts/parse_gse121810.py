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
Parse GSE121810 glioma immunotherapy dataset
"""
import os
import gzip
import re
import pandas as pd

OUT_DIR = os.path.join(ROOT, "output", "h2_bulk")

# Read expression matrix
print('Reading expression matrix...')
expr_file = os.path.join(OUT_DIR, 'GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx')
expr_df = pd.read_excel(expr_file, index_col=0)
print(f'Expression matrix shape: {expr_df.shape}')
print(f'Number of genes: {expr_df.shape[0]}')
print(f'Number of samples: {expr_df.shape[1]}')

# View sample names
print('\nSample name list:')
for i, col in enumerate(expr_df.columns):
    print(f'{i+1}: {col}')

# Parse sample information
print('\nParsing sample information...')
sample_info = []
for col in expr_df.columns:
    # Parse Pt3_A, Pt12_B, etc.
    parts = col.split('_')
    if len(parts) >= 2:
        patient_id = parts[0]
        group = parts[1]  # A or B
        sample_info.append({
            'sample_id': col,
            'patient_id': patient_id,
            'group': group
        })
    else:
        sample_info.append({
            'sample_id': col,
            'patient_id': col,
            'group': 'unknown'
        })

info_df = pd.DataFrame(sample_info)
print(f'\nSample information:')
print(info_df.to_string())

# Group statistics
print('\nGroup distribution:')
print(info_df['group'].value_counts())

# Check ZP3 expression
print('\nZP3 expression analysis:')
if 'ZP3' in expr_df.index:
    zp3_expr = expr_df.loc['ZP3']
    print(f'ZP3 expression values (first 10 samples):')
    for i, (sample, expr) in enumerate(zp3_expr.items()):
        if i < 10:
            print(f'{sample}: {expr:.3f}')
    
    # Statistics by group
    print('\nStatistics by group:')
    for group in info_df['group'].unique():
        samples = info_df[info_df['group'] == group]['sample_id'].values
        group_expr = zp3_expr[samples]
        print(f'Group {group} (n={len(samples)}): mean={group_expr.mean():.3f}, median={group_expr.median():.3f}, std={group_expr.std():.3f}')
else:
    print('ZP3 gene not found')
    # Search for possible ZP3-related genes
    zp3_candidates = [g for g in expr_df.index if 'ZP3' in str(g).upper()]
    if zp3_candidates:
        print(f'Found possible ZP3-related genes: {zp3_candidates}')

# Parse series matrix to get clinical information
print('\nParsing series matrix...')
series_file = os.path.join(OUT_DIR, 'GSE121810_series_matrix.txt.gz')
if os.path.exists(series_file):
    lines = []
    with gzip.open(series_file, 'rt') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!'):
                lines.append(line.rstrip())
    
    print(f'Read {len(lines)} lines of metadata')
    
    # Extract fields
    fields = {}
    for line in lines:
        m = re.match(r'!(\w+)\s+(.*)', line, re.S)
        if not m:
            continue
        key = m.group(1)
        cols = m.group(2).split('\t')
        cols = [c.strip('"') for c in cols]
        fields[key] = cols
    
    print(f'Extracted {len(fields)} fields')
    
    # Inspect key fields
    key_fields = ['Sample_description', 'Sample_characteristics_ch1', 'Sample_characteristics_ch2']
    for field in key_fields:
        if field in fields:
            values = fields[field]
            print(f'\n{field} (length {len(values)}):')
            for i, v in enumerate(values[:5]):
                print(f'  {i}: {v}')
            if len(values) > 5:
                print(f'  ... ({len(values)} total)')
    
    # Extract treatment group information
    if 'Sample_characteristics_ch1' in fields:
        print('\nTreatment group information:')
        for i, v in enumerate(fields['Sample_characteristics_ch1']):
            if i < 5:
                print(f'  {i}: {v}')
    
    # Extract survival information
    survival_keywords = ['survival', 'os', 'overall', 'time', 'month', 'day', 'year', 'status', 'event', 'alive', 'dead', 'death']
    for key, values in fields.items():
        if any(keyword in key.lower() for keyword in survival_keywords):
            print(f'\n{key}: {values[:5]}...')
else:
    print('series matrix file does not exist')

# Save sample info
info_csv = os.path.join(OUT_DIR, 'gse121810_sample_info.csv')
info_df.to_csv(info_csv, index=False)
print(f'\nSample info saved: {info_csv}')

# Save ZP3 expression
if 'ZP3' in expr_df.index:
    zp3_csv = os.path.join(OUT_DIR, 'gse121810_zp3_expression.csv')
    zp3_df = pd.DataFrame({'sample': zp3_expr.index, 'zp3_expression': zp3_expr.values})
    zp3_df.to_csv(zp3_csv, index=False)
    print(f'ZP3 expression saved: {zp3_csv}')

# Save expression matrix
expr_csv = os.path.join(OUT_DIR, 'gse121810_expression.csv')
expr_df.to_csv(expr_csv)
print(f'Expression matrix saved: {expr_csv}')
