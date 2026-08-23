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
解析 GSE121810 胶质瘤免疫治疗数据集
"""
import os
import gzip
import re
import pandas as pd

OUT_DIR = os.path.join(ROOT, "output", "h2_bulk")

# 读取表达矩阵
print('读取表达矩阵...')
expr_file = os.path.join(OUT_DIR, 'GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx')
expr_df = pd.read_excel(expr_file, index_col=0)
print(f'表达矩阵形状: {expr_df.shape}')
print(f'基因数: {expr_df.shape[0]}')
print(f'样本数: {expr_df.shape[1]}')

# 查看样本名
print('\n样本名列表:')
for i, col in enumerate(expr_df.columns):
    print(f'{i+1}: {col}')

# 解析样本信息
print('\n解析样本信息...')
sample_info = []
for col in expr_df.columns:
    # 解析 Pt3_A, Pt12_B 等
    parts = col.split('_')
    if len(parts) >= 2:
        patient_id = parts[0]
        group = parts[1]  # A 或 B
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
print(f'\n样本信息:')
print(info_df.to_string())

# 统计分组
print('\n分组分布:')
print(info_df['group'].value_counts())

# 查看 ZP3 表达
print('\nZP3 表达分析:')
if 'ZP3' in expr_df.index:
    zp3_expr = expr_df.loc['ZP3']
    print(f'ZP3 表达值 (前10个样本):')
    for i, (sample, expr) in enumerate(zp3_expr.items()):
        if i < 10:
            print(f'{sample}: {expr:.3f}')
    
    # 按分组统计
    print('\n按分组统计:')
    for group in info_df['group'].unique():
        samples = info_df[info_df['group'] == group]['sample_id'].values
        group_expr = zp3_expr[samples]
        print(f'组 {group} (n={len(samples)}): mean={group_expr.mean():.3f}, median={group_expr.median():.3f}, std={group_expr.std():.3f}')
else:
    print('未找到 ZP3 基因')
    # 查找可能的 ZP3 相关基因
    zp3_candidates = [g for g in expr_df.index if 'ZP3' in str(g).upper()]
    if zp3_candidates:
        print(f'找到可能的 ZP3 相关基因: {zp3_candidates}')

# 解析 series matrix 获取临床信息
print('\n解析 series matrix...')
series_file = os.path.join(OUT_DIR, 'GSE121810_series_matrix.txt.gz')
if os.path.exists(series_file):
    lines = []
    with gzip.open(series_file, 'rt') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!'):
                lines.append(line.rstrip())
    
    print(f'读取 {len(lines)} 行元数据')
    
    # 提取字段
    fields = {}
    for line in lines:
        m = re.match(r'!(\w+)\s+(.*)', line, re.S)
        if not m:
            continue
        key = m.group(1)
        cols = m.group(2).split('\t')
        cols = [c.strip('"') for c in cols]
        fields[key] = cols
    
    print(f'提取到 {len(fields)} 个字段')
    
    # 查看关键字段
    key_fields = ['Sample_description', 'Sample_characteristics_ch1', 'Sample_characteristics_ch2']
    for field in key_fields:
        if field in fields:
            values = fields[field]
            print(f'\n{field} (长度 {len(values)}):')
            for i, v in enumerate(values[:5]):
                print(f'  {i}: {v}')
            if len(values) > 5:
                print(f'  ... (共 {len(values)} 个)')
    
    # 提取治疗组信息
    if 'Sample_characteristics_ch1' in fields:
        print('\n治疗组信息:')
        for i, v in enumerate(fields['Sample_characteristics_ch1']):
            if i < 5:
                print(f'  {i}: {v}')
    
    # 提取生存信息
    survival_keywords = ['survival', 'os', 'overall', 'time', 'month', 'day', 'year', 'status', 'event', 'alive', 'dead', 'death']
    for key, values in fields.items():
        if any(keyword in key.lower() for keyword in survival_keywords):
            print(f'\n{key}: {values[:5]}...')
else:
    print('series matrix 文件不存在')

# 保存样本信息
info_csv = os.path.join(OUT_DIR, 'gse121810_sample_info.csv')
info_df.to_csv(info_csv, index=False)
print(f'\n样本信息已保存: {info_csv}')

# 保存 ZP3 表达
if 'ZP3' in expr_df.index:
    zp3_csv = os.path.join(OUT_DIR, 'gse121810_zp3_expression.csv')
    zp3_df = pd.DataFrame({'sample': zp3_expr.index, 'zp3_expression': zp3_expr.values})
    zp3_df.to_csv(zp3_csv, index=False)
    print(f'ZP3 表达已保存: {zp3_csv}')

# 保存表达矩阵
expr_csv = os.path.join(OUT_DIR, 'gse121810_expression.csv')
expr_df.to_csv(expr_csv)
print(f'表达矩阵已保存: {expr_csv}')