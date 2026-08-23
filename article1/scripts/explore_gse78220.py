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
探索 GSE78220 数据结构和 ZP3 表达
"""
import os
import pandas as pd
import numpy as np

OUT_DIR = os.path.join(ROOT, "article1", "results")
local_file = os.path.join(OUT_DIR, 'GSE78220_PatientFPKM.xlsx')

# 读取表达矩阵
print('读取表达矩阵...')
expr_df = pd.read_excel(local_file, sheet_name=0, index_col=0)
print(f'表达矩阵形状: {expr_df.shape}')
print(f'基因数: {expr_df.shape[0]}')
print(f'样本数: {expr_df.shape[1]}')

# 查看列名（样本名）
print('\n样本名列表:')
for i, col in enumerate(expr_df.columns):
    print(f'{i+1}: {col}')

# 解析样本信息
print('\n解析样本信息...')
sample_info = []
for col in expr_df.columns:
    # 解析 Pt1.baseline, Pt16.OnTx 等
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
print(f'\n样本信息:')
print(info_df.to_string())

# 统计时间点
print('\n时间点分布:')
print(info_df['timepoint'].value_counts())

# 查看 ZP3 表达
print('\nZP3 表达分析:')
if 'ZP3' in expr_df.index:
    zp3_expr = expr_df.loc['ZP3']
    print(f'ZP3 表达值 (前10个样本):')
    for i, (sample, expr) in enumerate(zp3_expr.items()):
        if i < 10:
            print(f'{sample}: {expr:.3f}')
    
    # 按时间点分组
    print('\n按时间点分组统计:')
    for tp in info_df['timepoint'].unique():
        samples = info_df[info_df['timepoint'] == tp]['sample_id'].values
        tp_expr = zp3_expr[samples]
        print(f'{tp} (n={len(samples)}): mean={tp_expr.mean():.3f}, median={tp_expr.median():.3f}, std={tp_expr.std():.3f}')
else:
    print('未找到 ZP3 基因')
    # 查找可能的 ZP3 相关基因
    zp3_candidates = [g for g in expr_df.index if 'ZP3' in str(g).upper()]
    if zp3_candidates:
        print(f'找到可能的 ZP3 相关基因: {zp3_candidates}')
    else:
        print('未找到 ZP3 相关基因')

# 查看其他可能相关的基因
print('\n其他免疫相关基因表达:')
immune_genes = ['TREM2', 'CD68', 'CD163', 'PD-L1', 'CD274', 'PD1', 'PDCD1', 'CTLA4', 'IDO1', 'ARG1']
for gene in immune_genes:
    if gene in expr_df.index:
        expr = expr_df.loc[gene]
        print(f'{gene}: mean={expr.mean():.3f}, median={expr.median():.3f}, range={expr.min():.3f}-{expr.max():.3f}')

# 保存样本信息
info_csv = os.path.join(OUT_DIR, 'gse78220_sample_info.csv')
info_df.to_csv(info_csv, index=False)
print(f'\n样本信息已保存: {info_csv}')

# 保存 ZP3 表达
zp3_csv = os.path.join(OUT_DIR, 'gse78220_zp3_expression.csv')
zp3_df = pd.DataFrame({'sample': zp3_expr.index, 'zp3_expression': zp3_expr.values})
zp3_df.to_csv(zp3_csv, index=False)
print(f'ZP3 表达已保存: {zp3_csv}')