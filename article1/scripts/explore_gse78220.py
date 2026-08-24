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
Explore GSE78220 data structure and ZP3 expression
"""
import os
import pandas as pd
import numpy as np

OUT_DIR = os.path.join(ROOT, "article1", "results")
local_file = os.path.join(OUT_DIR, 'GSE78220_PatientFPKM.xlsx')

# Read expression matrix
print('Reading expression matrix...')
expr_df = pd.read_excel(local_file, sheet_name=0, index_col=0)
print(f'Expression matrix shape: {expr_df.shape}')
print(f'Number of genes: {expr_df.shape[0]}')
print(f'Number of samples: {expr_df.shape[1]}')

# View column names (sample names)
print('\nSample name list:')
for i, col in enumerate(expr_df.columns):
    print(f'{i+1}: {col}')

# Parse sample information
print('\nParsing sample information...')
sample_info = []
for col in expr_df.columns:
    # Parse Pt1.baseline, Pt16.OnTx, etc.
    parts = col.split('.')
    if len(parts) >= 2:
        patient_id = parts[0]
        timepoint = parts[1]
        sample_info.append({
            'sample_id': col,
            'patient_id': patient_id,
            'timepoint': timepoint
        })
    else:
        sample_info.append({
            'sample_id': col,
            'patient_id': col,
            'timepoint': 'unknown'
        })

info_df = pd.DataFrame(sample_info)
print(f'\nSample information:')
print(info_df.to_string())

# Count timepoints
print('\nTimepoint distribution:')
print(info_df['timepoint'].value_counts())

# Check ZP3 expression
print('\nZP3 expression analysis:')
if 'ZP3' in expr_df.index:
    zp3_expr = expr_df.loc['ZP3']
    print(f'ZP3 expression values (first 10 samples):')
    for i, (sample, expr) in enumerate(zp3_expr.items()):
        if i < 10:
            print(f'{sample}: {expr:.3f}')
    
    # Group by timepoint
    print('\nStatistics grouped by timepoint:')
    for tp in info_df['timepoint'].unique():
        samples = info_df[info_df['timepoint'] == tp]['sample_id'].values
        tp_expr = zp3_expr[samples]
        print(f'{tp} (n={len(samples)}): mean={tp_expr.mean():.3f}, median={tp_expr.median():.3f}, std={tp_expr.std():.3f}')
else:
    print('ZP3 gene not found')
    # Search for possible ZP3-related genes
    zp3_candidates = [g for g in expr_df.index if 'ZP3' in str(g).upper()]
    if zp3_candidates:
        print(f'Found possible ZP3-related genes: {zp3_candidates}')
    else:
        print('No ZP3-related genes found')

# Check other potentially related genes
print('\nOther immune-related gene expression:')
immune_genes = ['TREM2', 'CD68', 'CD163', 'PD-L1', 'CD274', 'PD1', 'PDCD1', 'CTLA4', 'IDO1', 'ARG1']
for gene in immune_genes:
    if gene in expr_df.index:
        expr = expr_df.loc[gene]
        print(f'{gene}: mean={expr.mean():.3f}, median={expr.median():.3f}, range={expr.min():.3f}-{expr.max():.3f}')

# Save sample information
info_csv = os.path.join(OUT_DIR, 'gse78220_sample_info.csv')
info_df.to_csv(info_csv, index=False)
print(f'\nSample information saved: {info_csv}')

# Save ZP3 expression
zp3_csv = os.path.join(OUT_DIR, 'gse78220_zp3_expression.csv')
zp3_df = pd.DataFrame({'sample': zp3_expr.index, 'zp3_expression': zp3_expr.values})
zp3_df.to_csv(zp3_csv, index=False)
print(f'ZP3 expression saved: {zp3_csv}')
