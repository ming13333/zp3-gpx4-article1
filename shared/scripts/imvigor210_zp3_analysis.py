#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
IMvigor210 队列 ZP3 与免疫治疗反应分析
分析 ZP3 在尿路上皮癌免疫治疗中的预测价值
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import mannwhitneyu, chi2_contingency
import warnings
warnings.filterwarnings('ignore')

# 设置中文显示
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

print("=" * 60)
print("IMvigor210 队列 ZP3 免疫治疗反应分析")
print("=" * 60)

# ============================================================
# 1. 加载数据
# ============================================================
print("\n[1] 加载数据...")

# 加载表型数据
IMDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'output', 'immunotherapy_validation')
pData = pd.read_csv(os.path.join(IMDIR, 'pData_IMvigor210.csv'))
print(f"  表型数据: {pData.shape[0]} 样本, {pData.shape[1]} 变量")

# 加载特征数据（基因注释）
fData = pd.read_csv(os.path.join(IMDIR, 'fData_IMvigor210.csv'))
print(f"  特征数据: {fData.shape[0]} 基因")

# 加载表达矩阵
print("  加载表达矩阵（这可能需要一些时间）...")
exmat = pd.read_csv(os.path.join(IMDIR, 'exmat_censored_IMvigor210.csv'), index_col=0)
print(f"  表达矩阵: {exmat.shape[0]} 基因 × {exmat.shape[1]} 样本")

# ============ 显式行匹配校验（2026-08-17 科学审计修复）============
# exmat 行名为 gene_1..N（按 fData 顺序写入）。必须在分析前验证：
#  (a) 行数与 fData 一致；(b) ZP3 对应行（fData Symbol 列）能被定位且表达合理。
assert exmat.shape[0] == fData.shape[0], \
    f"行数不匹配: exmat {exmat.shape[0]} vs fData {fData.shape[0]}"
fData_tmp = fData.copy()
fData_tmp['row_name'] = ['gene_' + str(i) for i in range(1, len(fData_tmp) + 1)]
_sym2row = dict(zip(fData_tmp['Symbol'].astype(str), fData_tmp['row_name']))
_zp3_row = _sym2row.get('ZP3')
assert _zp3_row is not None, "fData 中未找到 ZP3 符号"
_zp3_val = exmat.loc[_zp3_row].astype(float)
assert (_zp3_val > 0).mean() > 0.5, f"ZP3 行非零比例过低: {(_zp3_val > 0).mean():.1%}"
print(f"  [行匹配校验通过] fData 行数={fData.shape[0]}, ZP3 位于行 '{_zp3_row}' (Entrez 7784), "
      f"非零比例={(_zp3_val > 0).mean():.1%}, 中位数={_zp3_val.median():.1f}")

# ============================================================
# 2. 数据预处理
# ============================================================
print("\n[2] 数据预处理...")

# 创建基因名映射（从fData）
gene_map = fData.set_index('entrez_id')['symbol'].to_dict()

# 将表达矩阵的行名从gene_XXX转换为实际基因名
# 表达矩阵的行名是gene_1, gene_2... 对应fData的顺序
exmat.index = [gene_map.get(fData.iloc[i, 0], f"gene_{i}") for i in range(len(exmat))]
print(f"  表达矩阵行名已更新为基因符号")

# 提取ZP3表达
if 'ZP3' in exmat.index:
    zp3_expr = exmat.loc['ZP3']
    print(f"  ZP3 表达已提取: {len(zp3_expr)} 样本")
    print(f"  ZP3 表达范围: {zp3_expr.min():.2f} - {zp3_expr.max():.2f}")
    print(f"  ZP3 中位表达: {zp3_expr.median():.2f}")
else:
    print("  错误: ZP3 基因未找到!")
    exit(1)

# 合并临床信息
pData_aligned = pData.set_index('X').loc[exmat.columns]
print(f"  对齐后样本数: {len(pData_aligned)}")

# ============================================================
# 3. ZP3 表达与治疗反应分析
# ============================================================
print("\n[3] ZP3 表达与治疗反应分析...")

# 添加ZP3表达到临床数据
pData_aligned['ZP3_expr'] = zp3_expr.values

# 定义响应组（排除NE - 不可评估）
responders = pData_aligned[pData_aligned['binaryResponse'] == 'CR/PR']
non_responders = pData_aligned[pData_aligned['binaryResponse'] == 'SD/PD']

print(f"  响应者 (CR/PR): {len(responders)} 例")
print(f"  非响应者 (SD/PD): {len(non_responders)} 例")

# ZP3 表达在响应组间的差异
zp3_resp = responders['ZP3_expr']
zp3_nonresp = non_responders['ZP3_expr']

# Mann-Whitney U 检验
stat_u, p_u = mannwhitneyu(zp3_resp, zp3_nonresp, alternative='two-sided')
print(f"\n  ZP3 表达差异检验 (Mann-Whitney U):")
print(f"    响应者 中位数: {zp3_resp.median():.2f} (IQR: {zp3_resp.quantile(0.25):.2f} - {zp3_resp.quantile(0.75):.2f})")
print(f"    非响应者 中位数: {zp3_nonresp.median():.2f} (IQR: {zp3_nonresp.quantile(0.25):.2f} - {zp3_nonresp.quantile(0.75):.2f})")
print(f"    U统计量: {stat_u:.2f}, p值: {p_u:.4f}")

# 效应量 (rank-biserial correlation)
# 2026-08-17 审计修复: 统一 r = 2U/(n1*n2) - 1 (U 为 responder 组), positive = responders 更高
n1, n2 = len(zp3_resp), len(zp3_nonresp)
r = (2 * stat_u) / (n1 * n2) - 1
print(f"    效应量 (rank-biserial r): {r:.4f}  (positive = higher ZP3 in responders)")

# ============================================================
# 4. ZP3 高/低表达分组分析
# ============================================================
print("\n[4] ZP3 高/低表达分组分析...")

# 使用中位数分割
median_zp3 = pData_aligned['ZP3_expr'].median()
pData_aligned['ZP3_group'] = pData_aligned['ZP3_expr'].apply(lambda x: 'High' if x >= median_zp3 else 'Low')

# 列联表
contingency = pd.crosstab(pData_aligned['ZP3_group'], pData_aligned['binaryResponse'])
print("\n  ZP3 分组 × 治疗反应 列联表:")
print(contingency)

# 卡方检验
chi2, p_chi, dof, expected = chi2_contingency(contingency)
print(f"\n  卡方检验: χ² = {chi2:.2f}, p = {p_chi:.4f}")

# 计算响应率
response_rates = pData_aligned.groupby('ZP3_group')['binaryResponse'].apply(
    lambda x: (x == 'CR/PR').sum() / len(x) * 100
)
print(f"\n  各组响应率:")
for group, rate in response_rates.items():
    print(f"    {group} ZP3: {rate:.1f}%")

# ============================================================
# 5. 生存分析
# ============================================================
print("\n[5] 生存分析...")

# 检查生存数据
if 'OS' in pData_aligned.columns and 'censored' in pData_aligned.columns:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    # KM曲线
    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()

    mask_high = pData_aligned['ZP3_group'] == 'High'
    mask_low = pData_aligned['ZP3_group'] == 'Low'

    kmf_high.fit(pData_aligned.loc[mask_high, 'OS'],
                 pData_aligned.loc[mask_high, 'censored'],
                 label='ZP3 High')
    kmf_low.fit(pData_aligned.loc[mask_low, 'OS'],
                pData_aligned.loc[mask_low, 'censored'],
                label='ZP3 Low')

    # Log-rank 检验
    results = logrank_test(pData_aligned.loc[mask_high, 'OS'],
                           pData_aligned.loc[mask_low, 'OS'],
                           pData_aligned.loc[mask_high, 'censored'],
                           pData_aligned.loc[mask_low, 'censored'])

    print(f"  中位生存期:")
    print(f"    ZP3 High: {kmf_high.median_survival_time_:.1f} 天")
    print(f"    ZP3 Low: {kmf_low.median_survival_time_:.1f} 天")
    print(f"  Log-rank 检验: χ² = {results.test_statistic:.2f}, p = {results.p_value:.4f}")

    # 绘制KM曲线
    fig, ax = plt.subplots(figsize=(10, 6))
    kmf_high.plot_survival_function(ax=ax, ci_show=True)
    kmf_low.plot_survival_function(ax=ax, ci_show=True)
    ax.set_title('IMvigor210: Overall Survival by ZP3 Expression', fontsize=14)
    ax.set_xlabel('Time (days)', fontsize=12)
    ax.set_ylabel('Survival Probability', fontsize=12)
    ax.legend(fontsize=12)
    ax.text(0.65, 0.95, f'Log-rank p = {results.p_value:.3f}',
            transform=ax.transAxes, fontsize=12,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    plt.tight_layout()
    plt.savefig('fig_imvigor210_km_zp3.png', dpi=300, bbox_inches='tight')
    print("  KM曲线已保存: fig_imvigor210_km_zp3.png")
else:
    print("  警告: 生存数据列不完整，跳过生存分析")

# ============================================================
# 6. 免疫特征相关性分析
# ============================================================
print("\n[6] 免疫特征相关性分析...")

# 定义免疫特征基因集
immune_signatures = {
    'T_cell_exhaustion': ['PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'TOX'],
    'Cytolytic_activity': ['GZMA', 'GZMB', 'PRF1', 'IFNG'],
    'Immune_checkpoint': ['CD274', 'PDCD1LG2', 'CTLA4', 'PDCD1', 'LAG3', 'HAVCR2'],
    'Macrophage_M2': ['CD163', 'MSR1', 'MRC1', 'VSIG4', 'CD206'],
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF18'],
    'TGF_beta': ['TGFB1', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2'],
    'IFN_gamma_response': ['STAT1', 'IRF1', 'CXCL10', 'CXCL9', 'IDO1']
}

# 计算每个特征的平均表达
print("\n  免疫特征与ZP3相关性:")
correlation_results = []

for sig_name, genes in immune_signatures.items():
    # 检查基因是否存在
    valid_genes = [g for g in genes if g in exmat.index]
    if len(valid_genes) >= 2:
        # 计算特征评分（平均表达）
        sig_score = exmat.loc[valid_genes].mean(axis=0)
        # 与ZP3计算Spearman相关
        rho, p_spear = stats.spearmanr(zp3_expr.values, sig_score.values)
        correlation_results.append({
            'Signature': sig_name,
            'Genes_found': len(valid_genes),
            'Spearman_rho': rho,
            'p_value': p_spear
        })
        sig_str = f"    {sig_name}: ρ = {rho:.3f}, p = {p_spear:.4f}"
        if p_spear < 0.05:
            sig_str += " *"
        print(sig_str)

# ============================================================
# 7. 可视化
# ============================================================
print("\n[7] 生成可视化...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 7.1 ZP3表达分布（按响应组）
ax1 = axes[0, 0]
data_plot = pData_aligned[pData_aligned['binaryResponse'].isin(['CR/PR', 'SD/PD'])]
sns.boxplot(data=data_plot, x='binaryResponse', y='ZP3_expr', ax=ax1,
            palette=['#2ecc71', '#e74c3c'])
sns.stripplot(data=data_plot, x='binaryResponse', y='ZP3_expr', ax=ax1,
              color='black', alpha=0.3, size=3)
ax1.set_title('ZP3 Expression by Treatment Response', fontsize=12)
ax1.set_xlabel('Response Group', fontsize=11)
ax1.set_ylabel('ZP3 Expression (log2)', fontsize=11)
ax1.text(0.5, 0.95, f'Mann-Whitney p = {p_u:.3f}',
         transform=ax1.transAxes, ha='center', va='top',
         fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# 7.2 ZP3高/低组响应率
ax2 = axes[0, 1]
response_rates_df = pd.DataFrame({
    'ZP3 Group': ['Low', 'High'],
    'Response Rate (%)': [response_rates.get('Low', 0), response_rates.get('High', 0)]
})
bars = ax2.bar(response_rates_df['ZP3 Group'], response_rates_df['Response Rate (%)'],
               color=['#3498db', '#e74c3c'], edgecolor='black')
ax2.set_title('Response Rate by ZP3 Expression', fontsize=12)
ax2.set_ylabel('Response Rate (CR/PR, %)', fontsize=11)
ax2.set_ylim(0, max(response_rates_df['Response Rate (%)']) * 1.3)
for bar, val in zip(bars, response_rates_df['Response Rate (%)']):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
             f'{val:.1f}%', ha='center', va='bottom', fontsize=10)
ax2.text(0.5, 0.95, f'χ² p = {p_chi:.3f}',
         transform=ax2.transAxes, ha='center', va='top',
         fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# 7.3 免疫特征相关性热图
ax3 = axes[1, 0]
if correlation_results:
    corr_df = pd.DataFrame(correlation_results)
    corr_df = corr_df.set_index('Signature')
    sns.heatmap(corr_df[['Spearman_rho']], annot=True, fmt='.3f', cmap='RdBu_r',
                center=0, ax=ax3, cbar_kws={'label': 'Spearman ρ'})
    ax3.set_title('ZP3 vs Immune Signatures', fontsize=12)
    ax3.set_ylabel('')

# 7.4 ZP3 表达分布（整体）
ax4 = axes[1, 1]
ax4.hist(pData_aligned['ZP3_expr'], bins=30, color='steelblue', edgecolor='black', alpha=0.7)
ax4.axvline(median_zp3, color='red', linestyle='--', linewidth=2, label=f'Median = {median_zp3:.2f}')
ax4.set_title('Distribution of ZP3 Expression', fontsize=12)
ax4.set_xlabel('ZP3 Expression (log2)', fontsize=11)
ax4.set_ylabel('Frequency', fontsize=11)
ax4.legend()

plt.suptitle('IMvigor210 Cohort: ZP3 and Immunotherapy Response', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('fig_imvigor210_zp3_analysis.png', dpi=300, bbox_inches='tight')
print("  主图已保存: fig_imvigor210_zp3_analysis.png")

# ============================================================
# 8. 保存结果
# ============================================================
print("\n[8] 保存结果...")

# 保存详细结果
results_summary = {
    'Cohort': 'IMvigor210',
    'Cancer_Type': 'Urothelial_Carcinoma',
    'Treatment': 'Atezolizumab (anti-PD-L1)',
    'Total_samples': len(pData_aligned),
    'Responders_CR_PR': len(responders),
    'NonResponders_SD_PD': len(non_responders),
    'ZP3_median_responders': zp3_resp.median(),
    'ZP3_median_nonresponders': zp3_nonresp.median(),
    'MannWhitney_U': stat_u,
    'MannWhitney_p': p_u,
    'Effect_size_r': r,
    'Chi2_p': p_chi,
    'Response_rate_Low_ZP3': response_rates.get('Low', 0),
    'Response_rate_High_ZP3': response_rates.get('High', 0)
}

# 保存为CSV
results_df = pd.DataFrame([results_summary])
results_df.to_csv(os.path.join(IMDIR, 'imvigor210_zp3_results.csv'), index=False)
print("  结果已保存: imvigor210_zp3_results.csv")

# 保存免疫特征相关性
if correlation_results:
    corr_df = pd.DataFrame(correlation_results)
    corr_df.to_csv(os.path.join(IMDIR, 'imvigor210_zp3_immune_correlations.csv'), index=False)
    print("  免疫相关性已保存: imvigor210_zp3_immune_correlations.csv")

# ============================================================
# 9. 总结
# ============================================================
print("\n" + "=" * 60)
print("分析总结")
print("=" * 60)

print(f"""
队列: IMvigor210 (尿路上皮癌, Atezolizumab 治疗)
样本量: {len(pData_aligned)} 例 (CR/PR: {len(responders)}, SD/PD: {len(non_responders)})

主要发现:
1. ZP3 表达差异:
   - 响应者 vs 非响应者: p = {p_u:.4f}
   - 效应量 r = {r:.4f}

2. 治疗反应率:
   - ZP3 Low: {response_rates.get('Low', 0):.1f}%
   - ZP3 High: {response_rates.get('High', 0):.1f}%
   - 卡方检验 p = {p_chi:.4f}

结论:
{"ZP3 高表达与较差的免疫治疗反应相关" if r > 0 and p_u < 0.05 else
 "ZP3 低表达与较差的免疫治疗反应相关" if r < 0 and p_u < 0.05 else
 "ZP3 表达与免疫治疗反应无显著关联"}
""")

print("\n分析完成!")