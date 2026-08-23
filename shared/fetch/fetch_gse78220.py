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
获取 GSE78220 胶质瘤免疫治疗数据集
数据来源：GEO FTP 补充文件 GSE78220_PatientFPKM.xlsx
"""
import os
import requests
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# 设置路径
OUT_DIR = os.path.join(ROOT, "output", "h2_bulk")
os.makedirs(OUT_DIR, exist_ok=True)

# 下载文件
url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx'
local_file = os.path.join(OUT_DIR, 'GSE78220_PatientFPKM.xlsx')

print('正在下载 GSE78220 数据...')
try:
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    with open(local_file, 'wb') as f:
        f.write(r.content)
    print(f'下载完成: {local_file} ({len(r.content)} bytes)')
except Exception as e:
    print(f'下载失败: {e}')
    exit(1)

# 读取 Excel 文件
print('\n正在读取 Excel 文件...')
try:
    # 尝试读取所有 sheet
    xl = pd.ExcelFile(local_file)
    print(f'Sheet 名称: {xl.sheet_names}')
    
    # 通常第一个 sheet 是表达矩阵，第二个是临床数据
    # 先查看每个 sheet 的前几行
    for sheet in xl.sheet_names:
        df = pd.read_excel(local_file, sheet_name=sheet, nrows=5)
        print(f'\n=== Sheet: {sheet} ===')
        print(f'形状: {df.shape}')
        print(f'列名: {list(df.columns)[:10]}...')
        print(df.head(2).to_string())
        
except Exception as e:
    print(f'读取失败: {e}')
    exit(1)

# 假设第一个 sheet 是表达矩阵，第二个是临床数据
# 需要根据实际结构调整
print('\n正在解析数据...')

# 读取完整数据
try:
    # 表达矩阵
    expr_df = pd.read_excel(local_file, sheet_name=0, index_col=0)
    print(f'表达矩阵形状: {expr_df.shape}')
    
    # 临床数据
    clin_df = pd.read_excel(local_file, sheet_name=1, index_col=0)
    print(f'临床数据形状: {clin_df.shape}')
    print(f'临床数据列: {list(clin_df.columns)}')
    
    # 保存为 CSV
    expr_csv = os.path.join(OUT_DIR, 'gse78220_expression.csv')
    clin_csv = os.path.join(OUT_DIR, 'gse78220_clinical.csv')
    
    expr_df.to_csv(expr_csv)
    clin_df.to_csv(clin_csv)
    
    print(f'\n已保存:')
    print(f'表达矩阵: {expr_csv}')
    print(f'临床数据: {clin_csv}')
    
    # 查看 ZP3 表达情况
    if 'ZP3' in expr_df.index:
        zp3_expr = expr_df.loc['ZP3']
        print(f'\nZP3 表达统计:')
        print(f'样本数: {len(zp3_expr)}')
        print(f'均值: {zp3_expr.mean():.3f}')
        print(f'中位数: {zp3_expr.median():.3f}')
        print(f'标准差: {zp3_expr.std():.3f}')
        print(f'范围: {zp3_expr.min():.3f} - {zp3_expr.max():.3f}')
    else:
        print('\n警告: 未找到 ZP3 基因')
        # 查找可能的 ZP3 相关基因
        zp3_candidates = [g for g in expr_df.index if 'ZP3' in str(g).upper()]
        if zp3_candidates:
            print(f'找到可能的 ZP3 相关基因: {zp3_candidates}')
    
    # 查看临床数据中的治疗反应信息
    print(f'\n临床数据预览:')
    print(clin_df.head(10).to_string())
    
    # 查找可能的治疗反应列
    response_cols = [c for c in clin_df.columns if any(keyword in str(c).lower() for keyword in ['response', 'outcome', 'treatment', 'anti', 'pd', 'pd-l1', 'immunotherapy'])]
    if response_cols:
        print(f'\n可能的治疗反应相关列: {response_cols}')
        for col in response_cols[:3]:  # 只显示前3个
            print(f'{col}: {clin_df[col].value_counts().to_dict()}')
    
    # 查找生存时间列
    survival_cols = [c for c in clin_df.columns if any(keyword in str(c).lower() for keyword in ['survival', 'os', 'overall', 'time', 'month', 'day', 'year'])]
    if survival_cols:
        print(f'\n可能的生存相关列: {survival_cols}')
    
    # 查找状态列（死亡/生存）
    status_cols = [c for c in clin_df.columns if any(keyword in str(c).lower() for keyword in ['status', 'event', 'alive', 'dead', 'death', 'vital'])]
    if status_cols:
        print(f'\n可能的状态相关列: {status_cols}')
        
except Exception as e:
    print(f'解析失败: {e}')
    import traceback
    traceback.print_exc()
    exit(1)

print('\n数据获取完成。')