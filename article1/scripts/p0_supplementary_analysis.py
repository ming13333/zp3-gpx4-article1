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
P0 补充分析:
1. ZP3 高/低分组的免疫特征差异 (Mann-Whitney U)
2. GSE121810 治疗组 (新辅助 vs 辅助) 的 ZP3 及免疫评分差异
3. 跨队列森林图: ZP3-免疫特征效应量
"""
import os, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

OUT = os.path.join(ROOT, "output", "h2_bulk")

# ============================================================
# 1. 加载已有数据
# ============================================================
scores78 = pd.read_csv(f'{OUT}/p0_gse78220_all_scores.csv', index_col=0)
scores12 = pd.read_csv(f'{OUT}/p0_gse121810_all_scores.csv', index_col=0)
corr_all = pd.read_csv(f'{OUT}/p0_zp3_correlation_all_features.csv')

print('GSE78220 scores:', scores78.shape)
print('GSE121810 scores:', scores12.shape)

# ============================================================
# 2. ZP3 高/低分组免疫特征差异
# ============================================================
print('\n' + '='*70)
print('  ZP3 HIGH vs LOW: 免疫特征差异')
print('='*70)

def zp3_group_analysis(df, cohort_name, pct=50):
    """ZP3 > median vs <= median 的免疫特征 Mann-Whitney U"""
    zp3 = df['ZP3'].values
    med = np.median(zp3)
    hi = zp3 > med
    lo = ~hi
    print(f'\n--- {cohort_name} (ZP3 median={med:.2f}, high={hi.sum()}, low={lo.sum()}) ---')

    results = []
    features = [c for c in df.columns if c not in ['ZP3']]
    for feat in features:
        v = df[feat].values
        mask = np.isfinite(v)
        if mask.sum() < 10: continue
        h = v[hi & mask]; l = v[lo & mask]
        if len(h) < 3 or len(l) < 3: continue
        stat, p = stats.mannwhitneyu(h, l, alternative='two-sided')
        fc = np.mean(h) / max(np.mean(l), 1e-10)
        results.append({
            'Feature': feat,
            'ZP3hi_mean': np.mean(h), 'ZP3lo_mean': np.mean(l),
            'log2FC': np.log2(fc), 'pvalue': p,
            'significant': p < 0.05
        })
    res_df = pd.DataFrame(results).sort_values('pvalue')
    # BH FDR
    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(res_df['pvalue'].values, method='fdr_bh')
    res_df['FDR'] = fdr
    res_df['BH_sig'] = fdr < 0.05

    # 打印显著结果
    sig = res_df[res_df['pvalue'] < 0.05]
    print(f'  p<0.05 特征: {len(sig)}/{len(res_df)}')
    if len(sig) > 0:
        for _, row in sig.iterrows():
            marker = '**' if row['BH_sig'] else '*'
            print(f'    {row["Feature"]:22s}: log2FC={row["log2FC"]:+.2f}  p={row["pvalue"]:.4f}  FDR={row["FDR"]:.4f} {marker}')

    return res_df

res78_grp = zp3_group_analysis(scores78, 'GSE78220 Melanoma')
res12_grp = zp3_group_analysis(scores12, 'GSE121810 GBM')

# ============================================================
# 3. GSE121810: 治疗组间免疫评分差异 (已有临床数据)
# ============================================================
print('\n' + '='*70)
print('  GSE121810: 治疗组 (新辅助 vs 辅助) 免疫评分差异')
print('='*70)

clinical12 = pd.read_csv(f'{OUT}/gse121810_sample_info.csv')
# 对齐样本
common = list(set(scores12.index) & set(clinical12.iloc[:,0].values))
if len(common) < 5:
    # 尝试从 annotation 匹配
    clinical12 = pd.read_csv(f'{OUT}/gse121810_all_fields.csv')
    print('  clinical columns:', list(clinical12.columns)[:10])
    # 尝试找 treatment group
    for c in clinical12.columns:
        if 'treat' in c.lower() or 'group' in c.lower() or 'type' in c.lower():
            print(f'  Found column: {c}')
            print(f'  Values: {clinical12[c].value_counts().to_dict()}')

# 从之前的分析我们知道新辅助 vs 辅助
# 重新读取临床数据
annot12 = pd.read_csv(f'{OUT}/gse121810_all_fields.csv')
print('  clinical12 columns:', list(annot12.columns))

# 检查是否有治疗组信息
therapy_col = None
for c in annot12.columns:
    vals = annot12[c].astype(str).str.lower()
    if any('neoadjuv' in v or 'adjuvant' in v or 'naive' in v for v in vals):
        therapy_col = c
        break

if therapy_col:
    print(f'\n  找到治疗组列: {therapy_col}')
    print(f'  分布: {annot12[therapy_col].value_counts().to_dict()}')
else:
    print('\n  未找到治疗组列，用样本名推断')
    # 从样本名推断: _A = adjuvant, _NA = neoadjuvant
    # 或者从之前的分析: Pt_A, Pt_NA 等
    sample_col = annot12.columns[0]
    print(f'  样本名列: {sample_col}')
    print(f'  样本: {annot12[sample_col].head().tolist()}')

# ============================================================
# 4. 可视化: ZP3 high vs low 对比图 + 森林图
# ============================================================
print('\n--- 生成可视化 ---')

fig = plt.figure(figsize=(18, 14))
gs = GridSpec(2, 2, hspace=0.35, wspace=0.3)

# --- 4a: GSE121810 ZP3 high vs low 免疫评分 (横向柱状图) ---
ax1 = fig.add_subplot(gs[0, 0])
# 选择有意义的特征
key_feats = ['Macrophages_M2', 'Macrophages_M1', 'Monocytes', 'CD163', 'MRC1',
             'Exhaustion', 'IFN_gamma', 'CYT', 'CD8_T', 'T_cells',
             'StromalScore', 'ImmuneScore', 'TREM2', 'IL10', 'TGFb', 'PD-L1']
key_feats = [f for f in key_feats if f in res12_grp['Feature'].values]
sub = res12_grp[res12_grp['Feature'].isin(key_feats)].copy()
sub = sub.sort_values('log2FC')

colors = ['#E24B4A' if v > 0 else '#378ADD' for v in sub['log2FC']]
bars = ax1.barh(range(len(sub)), sub['log2FC'].values, color=colors, edgecolor='white', linewidth=0.5)
ax1.set_yticks(range(len(sub)))
ax1.set_yticklabels(sub['Feature'].values, fontsize=10)
ax1.set_xlabel('log2FC (ZP3-high vs ZP3-low)', fontsize=11)
ax1.set_title('GSE121810 GBM: ZP3-high vs Low\nImmune Feature Differences', fontsize=13, fontweight='500')
ax1.axvline(0, color='gray', linestyle='--', linewidth=0.5)
# 标注显著性
for i, (_, row) in enumerate(sub.iterrows()):
    sig = '**' if row['BH_sig'] else '*' if row['pvalue'] < 0.05 else 'ns'
    ax1.text(row['log2FC'] + (0.02 if row['log2FC'] >= 0 else -0.02), i, sig,
             ha='left' if row['log2FC'] >= 0 else 'right', va='center', fontsize=9,
             fontweight='500', color='#791F1F' if sig != 'ns' else 'gray')

# --- 4b: GSE78220 ZP3 high vs low ---
ax2 = fig.add_subplot(gs[0, 1])
key_feats2 = [f for f in key_feats if f in res78_grp['Feature'].values]
sub2 = res78_grp[res78_grp['Feature'].isin(key_feats2)].copy()
sub2 = sub2.sort_values('log2FC')

colors2 = ['#E24B4A' if v > 0 else '#378ADD' for v in sub2['log2FC']]
ax2.barh(range(len(sub2)), sub2['log2FC'].values, color=colors2, edgecolor='white', linewidth=0.5)
ax2.set_yticks(range(len(sub2)))
ax2.set_yticklabels(sub2['Feature'].values, fontsize=10)
ax2.set_xlabel('log2FC (ZP3-high vs ZP3-low)', fontsize=11)
ax2.set_title('GSE78220 Melanoma: ZP3-high vs Low\nImmune Feature Differences', fontsize=13, fontweight='500')
ax2.axvline(0, color='gray', linestyle='--', linewidth=0.5)
for i, (_, row) in enumerate(sub2.iterrows()):
    sig = '**' if row['BH_sig'] else '*' if row['pvalue'] < 0.05 else 'ns'
    ax2.text(row['log2FC'] + (0.02 if row['log2FC'] >= 0 else -0.02), i, sig,
             ha='left' if row['log2FC'] >= 0 else 'right', va='center', fontsize=9,
             fontweight='500', color='#791F1F' if sig != 'ns' else 'gray')

# --- 4c: 森林图 (跨队列效应量) ---
ax3 = fig.add_subplot(gs[1, :])

# 合并两个队列的相关系数
merged = corr_all.merge(
    res78_grp[['Feature','log2FC']].rename(columns={'log2FC': 'log2FC_78'}),
    on='Feature', how='left'
).merge(
    res12_grp[['Feature','log2FC']].rename(columns={'log2FC': 'log2FC_12'}),
    on='Feature', how='left'
)

# 筛选有意义的特征 (|r|>0.3 在任一队列)
mask = (merged['Pearson_r_GSE78220'].abs() > 0.3) | (merged['Pearson_r_GSE121810'].abs() > 0.3)
forest = merged[mask].copy()
forest = forest.sort_values('Pearson_r_GSE121810', ascending=True)

y_pos = range(len(forest))
# GSE78220
ax3.errorbar(forest['Pearson_r_GSE78220'], y_pos, fmt='s', color='#378ADD',
             capsize=3, markersize=6, label='GSE78220 (Melanoma)', alpha=0.8)
# GSE121810
ax3.errorbar(forest['Pearson_r_GSE121810'], y_pos, fmt='o', color='#D85A30',
             capsize=3, markersize=6, label='GSE121810 (GBM)', alpha=0.8)

ax3.set_yticks(y_pos)
ax3.set_yticklabels(forest['Feature'].values, fontsize=10)
ax3.set_xlabel('Pearson r (ZP3 vs Feature)', fontsize=11)
ax3.set_title('Forest Plot: ZP3-Immune Feature Correlations\nCross-cohort Comparison', fontsize=13, fontweight='500')
ax3.axvline(0, color='gray', linestyle='--', linewidth=0.5)
ax3.axvline(0.3, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
ax3.axvline(-0.3, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
ax3.legend(loc='lower right', fontsize=10)
ax3.set_xlim(-0.75, 0.75)

plt.savefig(f'{OUT}/fig_p0_zp3_immunophenotype_analysis.png', dpi=200, bbox_inches='tight')
print(f'保存: fig_p0_zp3_immunophenotype_analysis.png')

# ============================================================
# 5. 保存所有结果
# ============================================================
res78_grp.to_csv(f'{OUT}/p0_gse78220_zp3_group_test.csv', index=False)
res12_grp.to_csv(f'{OUT}/p0_gse121810_zp3_group_test.csv', index=False)
print(f'保存: p0_gse78220_zp3_group_test.csv')
print(f'保存: p0_gse121810_zp3_group_test.csv')

print('\n=== P0 补充分析完成 ===')
