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
分析 GSE121810 中 ZP3 与治疗组的关系
数据：胶质瘤抗 PD-1 治疗队列（新辅助 vs 辅助 pembrolizumab）
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

OUT_DIR = os.path.join(ROOT, "article1", "results")

# 读取数据
print('读取数据...')
expr_df = pd.read_excel(os.path.join(OUT_DIR, 'GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx'), index_col=0)
annot_df = pd.read_csv(os.path.join(OUT_DIR, 'gse121810_sample_info.csv'))

# 提取治疗组信息
print('提取治疗组信息...')
# 从 series matrix 获取治疗组
series_file = os.path.join(OUT_DIR, 'gse121810_all_fields.csv')
series_df = pd.read_csv(series_file)

# 找到治疗信息列
therapy_col = None
for col in series_df.columns:
    if 'characteristics_ch1' in col:
        # 检查是否包含治疗信息
        values = series_df[col].dropna().astype(str)
        if any('therapy' in v.lower() for v in values):
            therapy_col = col
            break

if therapy_col:
    print(f'找到治疗信息列: {therapy_col}')
    therapy_info = series_df[therapy_col].apply(lambda x: 'neoadjuvant' if 'neoadjuvant' in str(x).lower() else 'adjuvant')
    print(f'治疗组分布: {therapy_info.value_counts().to_dict()}')
else:
    print('未找到治疗信息列，使用 A/B 分组')
    therapy_info = annot_df['group'].apply(lambda x: 'neoadjuvant' if x == 'A' else 'adjuvant')

# 合并数据
print('\n合并表达数据和临床信息...')
merged_data = []
for sample in expr_df.columns:
    if sample in annot_df['sample_id'].values:
        # 找到对应的样本信息
        sample_info = annot_df[annot_df['sample_id'] == sample].iloc[0]
        # 找到治疗组
        idx = annot_df[annot_df['sample_id'] == sample].index[0]
        therapy = therapy_info.iloc[idx] if idx < len(therapy_info) else 'unknown'
        
        merged_data.append({
            'sample_id': sample,
            'zp3_expr': expr_df.loc['ZP3', sample] if 'ZP3' in expr_df.index else np.nan,
            'trem2_expr': expr_df.loc['TREM2', sample] if 'TREM2' in expr_df.index else np.nan,
            'cd68_expr': expr_df.loc['CD68', sample] if 'CD68' in expr_df.index else np.nan,
            'cd163_expr': expr_df.loc['CD163', sample] if 'CD163' in expr_df.index else np.nan,
            'pdl1_expr': expr_df.loc['CD274', sample] if 'CD274' in expr_df.index else np.nan,
            'response': therapy,
            'group': sample_info['group'],
            'patient_id': sample_info['patient_id']
        })

merged_df = pd.DataFrame(merged_data)
print(f'合并后样本数: {len(merged_df)}')
print(f'治疗组分布: {merged_df["response"].value_counts().to_dict()}')

# 保存合并数据
merged_csv = os.path.join(OUT_DIR, 'gse121810_zp3_therapy.csv')
merged_df.to_csv(merged_csv, index=False)
print(f'合并数据已保存: {merged_csv}')

# 分析 ZP3 与治疗组的关系
print('\n=== ZP3 与治疗组的关系 ===')
print(f'\n按治疗组统计 ZP3 表达:')
therapy_stats = merged_df.groupby('response')['zp3_expr'].agg(['mean', 'median', 'std', 'count'])
print(therapy_stats)

# 统计检验
print('\n统计检验:')
neo_data = merged_df[merged_df['response'] == 'neoadjuvant']['zp3_expr'].dropna()
adj_data = merged_df[merged_df['response'] == 'adjuvant']['zp3_expr'].dropna()
if len(neo_data) > 0 and len(adj_data) > 0:
    t_stat, p_val = stats.ttest_ind(neo_data, adj_data)
    print(f'新辅助 vs 辅助: t={t_stat:.3f}, p={p_val:.4f}')
    
    # Mann-Whitney U 检验
    u_stat, p_mann = stats.mannwhitneyu(neo_data, adj_data)
    print(f'Mann-Whitney U: U={u_stat:.3f}, p={p_mann:.4f}')

# 可视化
print('\n生成可视化...')
plt.style.use('default')
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1. ZP3 表达按治疗组
ax1 = axes[0, 0]
therapy_order = ['neoadjuvant', 'adjuvant']
colors = {'neoadjuvant': '#3498db', 'adjuvant': '#e74c3c'}
for i, therapy in enumerate(therapy_order):
    data = merged_df[merged_df['response'] == therapy]['zp3_expr'].dropna()
    ax1.boxplot(data, positions=[i], widths=0.6, patch_artist=True,
                boxprops=dict(facecolor=colors[therapy], alpha=0.7),
                medianprops=dict(color='black', linewidth=2))
    # 添加散点
    jitter = np.random.default_rng(0).normal(0, 0.05, len(data))
    ax1.scatter([i] * len(data) + jitter, data, alpha=0.6, color=colors[therapy], s=30)
ax1.set_xticks(range(len(therapy_order)))
ax1.set_xticklabels(['Neoadjuvant', 'Adjuvant'])
ax1.set_ylabel('ZP3 Expression (FPKM)')
ax1.set_title('ZP3 Expression by Treatment Group')
ax1.grid(True, alpha=0.3)

# 添加统计检验结果
ax1.text(0.5, 0.95, f'p={p_val:.4f}', transform=ax1.transAxes, ha='center', va='top', fontsize=10)

# 2. ZP3 vs TREM2 散点图
ax2 = axes[0, 1]
for therapy in therapy_order:
    mask = merged_df['response'] == therapy
    ax2.scatter(merged_df[mask]['zp3_expr'], merged_df[mask]['trem2_expr'], 
                alpha=0.6, color=colors[therapy], label=therapy, s=50)
ax2.set_xlabel('ZP3 Expression')
ax2.set_ylabel('TREM2 Expression')
ax2.set_title('ZP3 vs TREM2 Expression')
ax2.legend()
ax2.grid(True, alpha=0.3)

# 3. ZP3 vs PD-L1 散点图
ax3 = axes[1, 0]
for therapy in therapy_order:
    mask = merged_df['response'] == therapy
    ax3.scatter(merged_df[mask]['zp3_expr'], merged_df[mask]['pdl1_expr'], 
                alpha=0.6, color=colors[therapy], label=therapy, s=50)
ax3.set_xlabel('ZP3 Expression')
ax3.set_ylabel('PD-L1 (CD274) Expression')
ax3.set_title('ZP3 vs PD-L1 Expression')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. ZP3 表达分布
ax4 = axes[1, 1]
neo_expr = merged_df[merged_df['response'] == 'neoadjuvant']['zp3_expr'].dropna()
adj_expr = merged_df[merged_df['response'] == 'adjuvant']['zp3_expr'].dropna()
ax4.hist(neo_expr, bins=10, alpha=0.6, color=colors['neoadjuvant'], label='Neoadjuvant', density=True)
ax4.hist(adj_expr, bins=10, alpha=0.6, color=colors['adjuvant'], label='Adjuvant', density=True)
ax4.set_xlabel('ZP3 Expression (FPKM)')
ax4.set_ylabel('Density')
ax4.set_title('ZP3 Expression Distribution')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
fig_path = os.path.join(OUT_DIR, 'fig_gse121810_zp3_therapy.png')
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.close()
print(f'图表已保存: {fig_path}')

# 计算相关系数
print('\n=== 相关性分析 ===')
corr_cols = ['zp3_expr', 'trem2_expr', 'cd68_expr', 'cd163_expr', 'pdl1_expr']
corr_df = merged_df[corr_cols].dropna()
if len(corr_df) > 5:
    corr_matrix = corr_df.corr()
    print('相关系数矩阵:')
    print(corr_matrix.round(3))
    
    # 保存相关矩阵
    corr_csv = os.path.join(OUT_DIR, 'gse121810_correlation_matrix.csv')
    corr_matrix.to_csv(corr_csv)
    print(f'相关矩阵已保存: {corr_csv}')

# 生成总结报告
print('\n=== 分析总结 ===')
summary = f"""
GSE121810 ZP3 与治疗组分析报告
==========================================

数据集: GSE121810 (胶质瘤抗 PD-1 治疗队列)
样本数: {len(merged_df)}
治疗组分布:
- 新辅助 pembrolizumab: {len(merged_df[merged_df['response'] == 'neoadjuvant'])} 例
- 辅助 pembrolizumab: {len(merged_df[merged_df['response'] == 'adjuvant'])} 例

ZP3 表达统计:
- 全部样本: mean={merged_df['zp3_expr'].mean():.3f}, median={merged_df['zp3_expr'].median():.3f}
- 新辅助组: mean={neo_data.mean():.3f}, median={neo_data.median():.3f}
- 辅助组: mean={adj_data.mean():.3f}, median={adj_data.median():.3f}

统计检验:
- 新辅助 vs 辅助: t={t_stat:.3f}, p={p_val:.4f}
- Mann-Whitney U: p={p_mann:.4f}

主要发现:
1. ZP3 在新辅助组和辅助组间的表达差异
2. ZP3 与 TREM2、PD-L1 的相关性
3. ZP3 作为潜在预测标志物的价值

图表: {fig_path}
"""
summary_path = os.path.join(OUT_DIR, 'gse121810_analysis_summary.txt')
with open(summary_path, 'w', encoding='utf-8') as f:
    f.write(summary)
print(f'总结报告已保存: {summary_path}')