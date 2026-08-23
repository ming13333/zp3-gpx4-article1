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
探索 GSE121810 的临床信息
"""
import os
import gzip
import re
import pandas as pd

OUT_DIR = os.path.join(ROOT, "article1", "results")
series_file = os.path.join(OUT_DIR, 'GSE121810_series_matrix.txt.gz')

print('解析 series matrix...')
lines = []
with gzip.open(series_file, 'rt') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            break
        if line.startswith('!'):
            lines.append(line.rstrip())

print(f'读取 {len(lines)} 行元数据')

# 提取所有字段
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

# 查看所有字段
print('\n所有字段:')
for key, values in fields.items():
    print(f'{key}: {len(values)} 个值')
    if len(values) <= 5:
        print(f'  值: {values}')
    else:
        print(f'  前5个值: {values[:5]}')

# 重点查看特征字段
print('\n=== 特征字段详情 ===')
characteristic_fields = [k for k in fields.keys() if k.startswith('Sample_characteristics')]
for field in characteristic_fields:
    values = fields[field]
    print(f'\n{field} (长度 {len(values)}):')
    # 统计唯一值
    unique_vals = list(set(values))
    print(f'  唯一值数量: {len(unique_vals)}')
    if len(unique_vals) <= 10:
        print(f'  唯一值: {unique_vals}')
    else:
        print(f'  示例: {unique_vals[:10]}')

# 提取治疗组信息
print('\n=== 治疗组信息 ===')
if 'Sample_characteristics_ch1' in fields:
    therapy_info = fields['Sample_characteristics_ch1']
    print(f'治疗信息: {therapy_info[:5]}...')
    
    # 统计治疗类型
    therapy_counts = pd.Series(therapy_info).value_counts()
    print(f'治疗类型分布:')
    for therapy, count in therapy_counts.items():
        print(f'  {therapy}: {count} 例')

# 提取生存信息
print('\n=== 生存信息 ===')
survival_keywords = ['survival', 'os', 'overall', 'time', 'month', 'day', 'year', 'status', 'event', 'alive', 'dead', 'death', 'progression', 'response', 'outcome']
for key, values in fields.items():
    if any(keyword in key.lower() for keyword in survival_keywords):
        print(f'{key}: {values[:5]}...')

# 提取患者ID映射
print('\n=== 患者ID映射 ===')
if 'Sample_title' in fields:
    titles = fields['Sample_title']
    descriptions = fields.get('Sample_description', [])
    print(f'样本标题 (前10个):')
    for i, title in enumerate(titles[:10]):
        desc = descriptions[i] if i < len(descriptions) else 'N/A'
        print(f'  {title}: {desc}')

# 保存所有字段
all_fields_csv = os.path.join(OUT_DIR, 'gse121810_all_fields.csv')
all_fields_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in fields.items()]))
all_fields_df.to_csv(all_fields_csv, index=False)
print(f'\n所有字段已保存: {all_fields_csv}')

# 尝试提取治疗反应信息
print('\n=== 尝试提取治疗反应信息 ===')
# 从 Sample_source_name_ch1 提取
if 'Sample_source_name_ch1' in fields:
    source_names = fields['Sample_source_name_ch1']
    print(f'样本来源 (前10个):')
    for i, source in enumerate(source_names[:10]):
        print(f'  {i}: {source}')
    
    # 统计来源类型
    source_counts = pd.Series(source_names).value_counts()
    print(f'来源类型分布:')
    for source, count in source_counts.items():
        print(f'  {source}: {count} 例')

# 从 Sample_characteristics_ch2 提取
if 'Sample_characteristics_ch2' in fields:
    ch2 = fields['Sample_characteristics_ch2']
    print(f'\n特征2 (前10个):')
    for i, val in enumerate(ch2[:10]):
        print(f'  {i}: {val}')
    
    # 统计唯一值
    unique_ch2 = list(set(ch2))
    print(f'特征2唯一值数量: {len(unique_ch2)}')
    if len(unique_ch2) <= 10:
        print(f'特征2唯一值: {unique_ch2}')