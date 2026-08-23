# -*- coding: utf-8 -*-
"""
P0-A1/A2/A3: 深度免疫表型分析 — GSE78220 + GSE121810
不依赖 xCell/MCPcounter/ESTIMATE 包，手工实现全部评分

产出:
  1. MCP-counter 风格细胞类型评分 (8种免疫 + 基质)
  2. ESTIMATE 风格评分 (immune/stromal/estimate)
  3. T cell 耗竭/CYT/IFN-γ 评分
  4. 完整免疫抑制因子面板
  5. ZP3 与所有评分的 Pearson/Spearman 相关矩阵
  6. 热图可视化
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

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output", "h2_bulk")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article1", "results")
np.random.seed(42)

# ============================================================
# 1. 基因集定义 (手工实现, 不依赖外部包)
# ============================================================

# --- MCP-counter 风格 (Becht 2016) ---
MCP_GENESETS = {
    'T_cells': ['CD2','CD3D','CD3E','CD3G','CD28','ICOS','LCK','ZAP70','TRAC','TRBC2'],
    'CD8_T': ['CD8A','CD8B','GZMA','GZMB','GZMK','PRF1','IFNG','NKG7','GNLY'],
    'NK_cells': ['NKG7','GNLY','KLRD1','KLRB1','KLRF1','NCAM1','KIR2DL1','KIR2DL3'],
    'B_lineage': ['CD79A','CD79B','MS4A1','CD19','CD20','CD180','TNFRSF17','SDC1'],
    'Monocytes': ['CD14','LYZ','S100A8','S100A9','FCN1','VCAN','CSF1R','CEBPA'],
    'Macrophages_M1': ['CD68','CD80','CD86','NOS2','IL1B','TNF','CXCL10','CD163'],
    'Macrophages_M2': ['CD163','MRC1','MSR1','VSIG4','FOLR2','TGFB1','IL10','ARG1','MARG1'],
    'Dendritic_cells': ['CD1C','CLEC9A','FCER1A','ITGAX','HLA-DRA','HLA-DQA1','LAMP3','BATF3'],
    'Endothelial': ['PECAM1','VWF','CDH5','ERG','KDR','FLT1','NRP1','ANGPT2'],
    'Fibroblasts': ['COL1A1','COL1A2','COL3A1','FAP','DCN','LUM','ACTA2','PDGFRA','PDGFRB'],
}

# --- ESTIMATE 风格 (Yoshihara 2013) ---
ESTIMATE_STROMAL = ['DCN','COL1A1','COL1A2','COL3A1','COL5A1','COL6A1','FAP','LUM',
                    'POSTN','FBLN1','FBLN2','FBLN5','BGN','OGN','ASPN','COL11A1',
                    'COL11A2','COL5A2','COL5A3','ADAM12','SULF1','SFRP1','SFRP2',
                    'WNT2','WIF1','DKK1','GREM1','CTHRC1','MFAP5','COMP','THBS2',
                    'SPARC','SPARCL1','SOD3','ISLR','ISLR2','ADH1B','ABCA8',
                    'PDGFRB','PDGFRA','RGS5','CSPG4','NOTCH3','TCF21','PDGFD']

ESTIMATE_IMMUNE = ['BIRC3','C1QA','C1QB','C1QC','CCL2','CCL3','CCL4','CCL5','CCL8',
                   'CCR5','CD14','CD163','CD2','CD3D','CD3E','CD3G','CD48','CD52',
                   'CD53','CD68','CD69','CD74','CD80','CD86','CD96','CLEC10A','CTSS',
                   'CXCL10','CXCL11','CXCL9','FCER1A','FCGR1A','GZMA','GZMB','GZMK',
                   'GZMN','HCK','ICOS','IDO1','IFNG','IL10','IL1B','IL2RG','IL6',
                   'IRF1','ITGAM','ITGAX','LCK','LY86','LYZ','NKG7','PDCD1','PTPRC',
                   'SLAMF7','SPN','TLR8','TNF','TXK']

# --- T cell 功能/耗竭 ---
EXHAUSTION_GENES = ['PDCD1','CTLA4','HAVCR2','LAG3','TIGIT','LAG3','TIGIT',
                     'TOX','TOX2','ENTPD1','LAYN','MAF','TGFB1','LAG3']
EXHAUSTION_GENES = list(dict.fromkeys(EXHAUSTION_GENES))  # deduplicate

CYT_GENES = ['GZMA','PRF1']

IFNG_SIGNATURE = ['IFNG','STAT1','CXCL9','CXCL10','CXCL11','IDO1','GBP1','GBP2',
                  'IRF1','HLA-DRA','HLA-DRB1','PSMB8','PSMB9','TAP1','TAP2']

# --- 完整免疫抑制因子 ---
IMMUNOSUPPRESSIVE_PANEL = {
    # 检查点
    'PD-L1': ['CD274','PDCD1LG2'],
    'CTLA4': ['CTLA4'],
    'TIM3': ['HAVCR2'],
    'LAG3': ['LAG3'],
    'TIGIT': ['TIGIT'],
    'VISTA': ['VSIR'],
    'B7-H3': ['CD276'],
    'BTLA': ['BTLA'],
    # 抑制性细胞因子
    'TGFb': ['TGFB1','TGFB2','TGFB3'],
    'IL10': ['IL10','IL10RA','IL10RB'],
    'VEGF': ['VEGFA','VEGFB','VEGFC','KDR','FLT1'],
    # 抑制性酶
    'IDO1': ['IDO1','IDO2'],
    'ARG1': ['ARG1','ARG2'],
    # 趋化因子 (招募免疫抑制细胞)
    'CCL2': ['CCL2'],
    'CCL5': ['CCL5'],
    'CXCL8': ['CXCL8','CXCL1'],
    # M2 标志物
    'CD163': ['CD163'],
    'MRC1': ['MRC1'],
    'MSR1': ['MSR1'],
    'VSIG4': ['VSIG4'],
    'FOLR2': ['FOLR2'],
    'TREM2': ['TREM2'],
    'MERTK': ['MERTK'],
    # 代谢抑制
    'CD39': ['ENTPD1'],
    'CD73': ['NT5E'],
    'B2M': ['B2M'],
}

def score_geneset(expr_df, gene_list, method='zscore'):
    """对表达矩阵计算基因集评分 (仅用存在的基因)。

    说明（回应实证审查「免疫去卷积为基因集均值简化版」）：
    本实现是【基因集 signature score】，并非绝对细胞比例反卷积
    （CIBERSORT/MCP-counter 本身也是基于参考谱的估计，非金标准）。
    两种评分模式：
      - 'zscore'（默认，推荐）：每基因先跨样本 z 标准化，再取基因集均值。
        优点：消除基因间量纲/基线差异，使不同基因集、不同队列间分数可比；
        高值代表该样本的基因集相对高表达（富集），低值代表相对低表达。
      - 'mean'：原始基因集表达均值（未标准化，仅供对照）。
    """
    avail = [g for g in gene_list if g in expr_df.index]
    if len(avail) == 0:
        return pd.Series(np.nan, index=expr_df.columns)
    sub = expr_df.loc[avail]
    if method == 'zscore':
        mu = sub.mean(axis=1)
        sd = sub.std(axis=1).replace(0, np.nan)
        z = sub.sub(mu, axis=0).div(sd, axis=0)
        z = z.replace([np.inf, -np.inf], np.nan)
        return z.mean(axis=0)
    if method == 'median':
        return sub.median(axis=0)
    return sub.mean(axis=0)  # 'mean'

def analyze_cohort(expr_path, clinical_path, cohort_name, expr_kwargs=None):
    """对单个队列做完整免疫表型分析"""
    print(f'\n{"="*70}')
    print(f'  {cohort_name}')
    print(f'{"="*70}')

    # 读表达矩阵
    if expr_kwargs is None:
        expr_kwargs = {}
    if expr_path.endswith('.xlsx'):
        expr = pd.read_excel(expr_path, index_col=0, **expr_kwargs)
    else:
        expr = pd.read_csv(expr_path, index_col=0, **expr_kwargs)
    # 转为 float
    expr = expr.apply(pd.to_numeric, errors='coerce').fillna(0)
    print(f'表达矩阵: {expr.shape[0]} 基因 x {expr.shape[1]} 样本')

    # 读临床数据
    clinical = pd.read_csv(clinical_path)

    # ---- A1: MCP-counter 风格 ----
    print('\n--- A1: MCP-counter 风格细胞评分 ---')
    mcp_scores = {}
    for name, genes in MCP_GENESETS.items():
        score = score_geneset(expr, genes)
        mcp_scores[name] = score
        avail = [g for g in genes if g in expr.index]
        print(f'  {name:22s}: avail={len(avail):2d}/{len(genes):2d}, mean={score.mean():.3f}, std={score.std():.3f}')
    mcp_df = pd.DataFrame(mcp_scores)

    # ---- A1b: ESTIMATE 风格 ----
    print('\n--- A1b: ESTIMATE 风格评分 ---')
    est_scores = {}
    for name, genes in [('StromalScore', ESTIMATE_STROMAL), ('ImmuneScore', ESTIMATE_IMMUNE)]:
        score = score_geneset(expr, genes)
        est_scores[name] = score
        avail = [g for g in genes if g in expr.index]
        print(f'  {name:20s}: avail={len(avail):2d}/{len(genes):2d}, mean={score.mean():.3f}, std={score.std():.3f}')
    est_scores['ESTIMATEScore'] = est_scores['StromalScore'] + est_scores['ImmuneScore']
    est_df = pd.DataFrame(est_scores)

    # ---- A2: T cell 功能/耗竭评分 ----
    print('\n--- A2: T cell 功能/耗竭评分 ---')
    func_scores = {}
    for name, genes in [('Exhaustion', EXHAUSTION_GENES), ('CYT', CYT_GENES), ('IFN_gamma', IFNG_SIGNATURE)]:
        score = score_geneset(expr, genes)
        func_scores[name] = score
        avail = [g for g in genes if g in expr.index]
        print(f'  {name:16s}: avail={len(avail):2d}/{len(genes):2d}, mean={score.mean():.3f}, std={score.std():.3f}')
    func_df = pd.DataFrame(func_scores)

    # ---- A3: 免疫抑制因子面板 ----
    print('\n--- A3: 免疫抑制因子面板 ---')
    immuno_scores = {}
    for name, genes in IMMUNOSUPPRESSIVE_PANEL.items():
        score = score_geneset(expr, genes)
        immuno_scores[name] = score
        avail = [g for g in genes if g in expr.index]
    immuno_df = pd.DataFrame(immuno_scores)
    # 打印摘要
    for name in immuno_df.columns:
        v = immuno_df[name]
        print(f'  {name:12s}: mean={v.mean():.3f}, std={v.std():.3f}')

    # ---- ZP3 表达 ----
    zp3 = score_geneset(expr, ['ZP3'])
    trem2 = score_geneset(expr, ['TREM2'])
    print(f'\nZP3: mean={zp3.mean():.3f}, std={zp3.std():.3f}')
    print(f'TREM2: mean={trem2.mean():.3f}, std={trem2.std():.3f}')

    # ---- 相关矩阵 ----
    print('\n--- ZP3 相关矩阵 ---')
    all_scores = pd.concat([mcp_df, est_df, func_df, immuno_df], axis=1)
    all_scores['ZP3'] = zp3
    all_scores['TREM2'] = trem2

    # Pearson + Spearman
    pearson_r = []
    pearson_p = []
    spearman_r = []
    spearman_p = []
    features = list(all_scores.columns)

    for feat in features:
        if feat == 'ZP3':
            pearson_r.append(1.0); pearson_p.append(0.0)
            spearman_r.append(1.0); spearman_p.append(0.0)
            continue
        x = all_scores['ZP3'].values
        y = all_scores[feat].values
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 5:
            pearson_r.append(np.nan); pearson_p.append(np.nan)
            spearman_r.append(np.nan); spearman_p.append(np.nan)
            continue
        pr, pp = stats.pearsonr(x[mask], y[mask])
        sr, sp = stats.spearmanr(x[mask], y[mask])
        pearson_r.append(pr); pearson_p.append(pp)
        spearman_r.append(sr); spearman_p.append(sp)

    corr_df = pd.DataFrame({
        'Feature': features,
        'Pearson_r': pearson_r,
        'Pearson_p': pearson_p,
        'Spearman_r': spearman_r,
        'Spearman_p': spearman_p,
    })

    # 打印显著相关
    print('\n  ZP3 vs 各特征 (Pearson, |r|>0.3 或 p<0.05):')
    for _, row in corr_df.iterrows():
        f = row['Feature']
        if f == 'ZP3': continue
        r = row['Pearson_r']
        p = row['Pearson_p']
        sig = '' if p > 0.05 else '*' if p > 0.01 else '**' if p > 0.001 else '***'
        if abs(r) > 0.3 or p < 0.05:
            print(f'    {f:22s}: r={r:+.3f} p={p:.4f} {sig}')

    return {
        'expr': expr, 'mcp_df': mcp_df, 'est_df': est_df,
        'func_df': func_df, 'immuno_df': immuno_df,
        'zp3': zp3, 'trem2': trem2, 'corr_df': corr_df,
        'all_scores': all_scores,
    }

# ============================================================
# 2. 运行分析
# ============================================================

# GSE78220 (黑色素瘤, n=28)
res78 = analyze_cohort(
    expr_path=f'{OUT}/GSE78220_PatientFPKM.xlsx',
    clinical_path=f'{OUT}/gse78220_annotation.csv',
    cohort_name='GSE78220 Melanoma anti-PD-1 (n=28)',
    expr_kwargs={'sheet_name': 0}
)

# GSE121810 (胶质瘤, n=29)
res12 = analyze_cohort(
    expr_path=f'{OUT}/gse121810_expression.csv',
    clinical_path=f'{OUT}/gse121810_sample_info.csv',
    cohort_name='GSE121810 GBM anti-PD-1 (n=29)'
)

# ============================================================
# 3. 跨队列合并相关矩阵 + 保存
# ============================================================
print('\n' + '='*70)
print('  跨队列 ZP3 相关对比')
print('='*70)

merged = res78['corr_df'][['Feature','Pearson_r','Pearson_p']].merge(
    res12['corr_df'][['Feature','Pearson_r','Pearson_p']],
    on='Feature', suffixes=('_GSE78220','_GSE121810')
)

# 过滤有意义的行 (|r|>0.2 在任一队列)
mask = (merged['Pearson_r_GSE78220'].abs() > 0.2) | (merged['Pearson_r_GSE121810'].abs() > 0.2)
merged_filtered = merged[mask].sort_values('Pearson_r_GSE78220', ascending=False)
print(merged_filtered.to_string(index=False, float_format='%.3f'))

merged.to_csv(f'{OUT}/p0_zp3_correlation_all_features.csv', index=False)
merged_filtered.to_csv(f'{OUT}/p0_zp3_correlation_filtered.csv', index=False)
print(f'\n保存: p0_zp3_correlation_all_features.csv')
print(f'保存: p0_zp3_correlation_filtered.csv')

# ============================================================
# 4. 可视化: 双队列相关热图
# ============================================================
print('\n--- 生成可视化 ---')

# 筛选有意义特征
keep = merged_filtered['Feature'].tolist()
if 'ZP3' not in keep: keep.append('ZP3')
if 'TREM2' not in keep: keep.append('TREM2')

# 构建热图矩阵
filtered_idx = merged_filtered.set_index('Feature')
# 确保 keep 中的 feature 都在 index 中
keep_valid = [f for f in keep if f in filtered_idx.index]
heat_data = pd.DataFrame({
    'GSE78220\n(Melanoma)': filtered_idx.loc[keep_valid, 'Pearson_r_GSE78220'],
    'GSE121810\n(GBM)': filtered_idx.loc[keep_valid, 'Pearson_r_GSE121810'],
}).T

fig, axes = plt.subplots(1, 2, figsize=(16, 10), gridspec_kw={'width_ratios': [1, 1.1]})

# 左: 热图
ax = axes[0]
cmap = plt.cm.RdBu_r
norm = mcolors.TwoSlopeNorm(vmin=-0.8, vcenter=0, vmax=0.8)
im = ax.imshow(heat_data.values, cmap=cmap, norm=norm, aspect='auto')
ax.set_xticks(range(len(heat_data.columns)))
ax.set_xticklabels(heat_data.columns, rotation=45, ha='right', fontsize=9)
ax.set_yticks(range(len(heat_data.index)))
ax.set_yticklabels(heat_data.index, fontsize=10)
ax.set_title('ZP3 vs Immune Features\nPearson r (cross-cohort)', fontsize=13, fontweight='500')
plt.colorbar(im, ax=ax, shrink=0.8, label='Pearson r')

# 在格子内标注数值
for i in range(len(heat_data.index)):
    for j in range(len(heat_data.columns)):
        val = heat_data.values[i, j]
        if not np.isnan(val):
            color = 'white' if abs(val) > 0.5 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=8, color=color)

# 右: 散点图 (GSE78220 r vs GSE121810 r)
ax2 = axes[1]
r78 = merged_filtered['Pearson_r_GSE78220'].values
r12 = merged_filtered['Pearson_r_GSE121810'].values
labels = merged_filtered['Feature'].values

ax2.scatter(r78, r12, c='steelblue', s=40, alpha=0.7, edgecolors='white', linewidth=0.5)
ax2.axhline(0, color='gray', linestyle='--', linewidth=0.5)
ax2.axvline(0, color='gray', linestyle='--', linewidth=0.5)
ax2.plot([-1, 1], [-1, 1], color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
ax2.set_xlabel('Pearson r in GSE78220 (Melanoma)', fontsize=11)
ax2.set_ylabel('Pearson r in GSE121810 (GBM)', fontsize=11)
ax2.set_title('Cross-cohort Reproducibility\nof ZP3-Immune Correlations', fontsize=13, fontweight='500')
ax2.set_xlim(-0.9, 0.9)
ax2.set_ylim(-0.9, 0.9)

    # 标注 r > 0.3 的点
for i, (x, y, lab) in enumerate(zip(r78, r12, labels)):
    if np.isfinite(x) and np.isfinite(y) and (abs(x) > 0.35 or abs(y) > 0.35):
        ax2.annotate(lab, (x, y), fontsize=7.5, ha='center', va='bottom',
                     xytext=(0, 4), textcoords='offset points')

# 计算跨队列相关
valid = np.isfinite(r78) & np.isfinite(r12)
cross_r, cross_p = stats.pearsonr(r78[valid], r12[valid])
ax2.text(0.05, 0.95, f'Cross-cohort r={cross_r:.3f}\np={cross_p:.4f}',
         transform=ax2.transAxes, fontsize=10, va='top',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig(f'{OUT}/fig_p0_zp3_immunophenotype_heatmap.png', dpi=200, bbox_inches='tight')
print(f'保存: fig_p0_zp3_immunophenotype_heatmap.png')

# ============================================================
# 5. 保存每个队列的完整评分表
# ============================================================
for name, res in [('gse78220', res78), ('gse121810', res12)]:
    res['all_scores'].to_csv(f'{OUT}/p0_{name}_all_scores.csv')
    res['corr_df'].to_csv(f'{OUT}/p0_{name}_zp3_correlation_full.csv', index=False)
    print(f'保存: p0_{name}_all_scores.csv, p0_{name}_zp3_correlation_full.csv')

print('\n=== P0 完成 ===')
