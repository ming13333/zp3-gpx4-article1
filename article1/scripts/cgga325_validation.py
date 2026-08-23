#!/usr/bin/env python3
"""
CGGA-325 胶质瘤独立验证集 ZP3 分析
与CGGA-693使用完全一致的分析逻辑和免疫基因集，
用于验证CGGA-693的发现是否可重复。
同时输出两队列的对比结果。
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import warnings
import os

warnings.filterwarnings('ignore')
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

base_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'output', 'cgga_validation')
output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'article1', 'results')
os.makedirs(output_dir, exist_ok=True)

# ========== 免疫基因集（与CGGA-693完全一致） ==========
immune_gene_sets = {
    'M2_Macrophage': ['CD163', 'MSR1', 'MRC1', 'CD206', 'TGFB1', 'IL10', 'ARG1', 'IDO1', 'CCL2', 'CCL5'],
    'T_cell_exhaustion': ['PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'BTLA', 'VSIR', 'IDO1', 'ENTPD1'],
    'Cytolytic_activity': ['PRF1', 'GZMA', 'GZMB', 'GZMH', 'GZMM', 'GNLY', 'NKG7', 'KLRK1', 'KLRD1', 'KLRC1'],
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF18', 'ICOS', 'CD274', 'PDCD1LG2', 'TIGIT', 'IL10', 'TGFB1'],
    'IFN_gamma': ['IFNG', 'STAT1', 'IRF1', 'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'GBP1', 'GBP2', 'PSMB8'],
    'Checkpoint': ['CD274', 'PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'BTLA', 'VSIR', 'IDO1', 'ENTPD1'],
    'Myeloid': ['CD33', 'CD14', 'ITGAM', 'CSF1R', 'MPO', 'ELANE', 'AZU1', 'PRTN3', 'CEACAM8', 'CEACAM6']
}

def load_and_analyze(clinical_path, expr_path, cohort_name):
    """加载并分析单个CGGA队列"""
    print(f"\n{'='*60}")
    print(f"  {cohort_name} 队列分析")
    print(f"{'='*60}")
    
    clinical = pd.read_csv(clinical_path, sep='\t')
    expr = pd.read_csv(expr_path, sep='\t', index_col=0)
    
    print(f"临床样本数: {clinical.shape[0]}")
    print(f"基因数: {expr.shape[0]}")
    
    zp3_expr = expr.loc['ZP3']
    common = sorted(set(clinical['CGGA_ID']) & set(zp3_expr.index))
    clinical_aligned = clinical[clinical['CGGA_ID'].isin(common)].set_index('CGGA_ID').loc[common]
    zp3_aligned = zp3_expr[common]
    
    print(f"有效样本数: {len(common)}")
    print(f"ZP3检出率: {(zp3_aligned>0).sum()}/{len(zp3_aligned)}")
    
    results = {'cohort': cohort_name, 'n': len(common), 'zp3': zp3_aligned, 'clinical': clinical_aligned}
    
    # ---- 免疫特征评分 ----
    immune_scores = {}
    for set_name, genes in immune_gene_sets.items():
        avail = [g for g in genes if g in expr.index]
        if len(avail) >= 3:
            scores = expr.loc[avail, common].mean(axis=0)
            immune_scores[set_name] = scores
    
    results['immune_scores'] = immune_scores
    return results

def compare_queue(r1, r2):
    """对比两个队列的ZP3-免疫特征相关性"""
    print(f"\n{'='*60}")
    print("  CGGA-693 vs CGGA-325 免疫特征相关性对比")
    print(f"{'='*60}")
    
    rows = []
    for name in r1['immune_scores'].keys():
        if name in r2['immune_scores']:
            rho1, p1 = stats.spearmanr(r1['zp3'], r1['immune_scores'][name][r1['clinical'].index])
            rho2, p2 = stats.spearmanr(r2['zp3'], r2['immune_scores'][name][r2['clinical'].index])
            rows.append({
                'Immune_feature': name,
                'CGGA693_rho': rho1, 'CGGA693_p': p1,
                'CGGA325_rho': rho2, 'CGGA325_p': p2,
            })
    
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(output_dir, 'cgga_693_vs_325_immune.csv'), index=False)
    
    # 打印表格
    print(f"\n{'特征':<18}{'693_rho':>8}{'693_p':>12}{'325_rho':>8}{'325_p':>12}")
    print('-'*60)
    for _, row in df.iterrows():
        sig1 = '*' if row['CGGA693_p'] < 0.05 else ' '
        sig2 = '*' if row['CGGA325_p'] < 0.05 else ' '
        print(f"{row['Immune_feature']:<16}{row['CGGA693_rho']:>8.3f}{row['CGGA693_p']:>12.2e}{row['CGGA325_rho']:>8.3f}{row['CGGA325_p']:>12.2e}")
    return df

def main():
    print("进行 CGGA-325 独立验证...")
    
    # 分析两个队列
    r693 = load_and_analyze(
        os.path.join(base_dir, 'CGGA.mRNAseq_693_clinical.20200506.txt'),
        os.path.join(base_dir, 'CGGA.mRNAseq_693.RSEM-genes.20200506.txt'),
        'CGGA-693'
    )
    r325 = load_and_analyze(
        os.path.join(base_dir, 'CGGA.mRNAseq_325_clinical.20200506.txt'),
        os.path.join(base_dir, 'CGGA.mRNAseq_325.RSEM-genes.20200506.txt'),
        'CGGA-325'
    )
    
    # 对比免疫特征
    comp = compare_queue(r693, r325)
    
    # CGGA-325 临床关联检验
    print(f"\n{'='*60}")
    print("  CGGA-325 临床特征关联")
    print(f"{'='*60}")
    
    clin = r325['clinical']
    zp3 = r325['zp3']
    
    # IDH
    idh_wt = clin[clin['IDH_mutation_status']=='Wildtype'].index
    idh_mut = clin[clin['IDH_mutation_status']=='Mutant'].index
    if len(idh_wt)>=5 and len(idh_mut)>=5:
        u, p_idh = stats.mannwhitneyu(zp3[idh_wt], zp3[idh_mut])
        print(f"IDH 野生型 vs 突变型: p={p_idh:.4f}")
        print(f"  IDH-WT 中位: {zp3[idh_wt].median():.3f}, IDH-Mut 中位: {zp3[idh_mut].median():.3f}")
    
    # Grade
    grade_grp = []
    grade_lab = []
    for g in ['WHO II','WHO III','WHO IV']:
        s = clin[clin['Grade']==g].index
        if len(s)>=10:
            grade_grp.append(zp3[s].values)
            grade_lab.append(g)
    if len(grade_grp)>=2:
        h, p_grade = stats.kruskal(*grade_grp)
        print(f"WHO分级: Kruskal-Wallis p={p_grade:.4f}")
    
    # 1p19q
    codel = clin[clin['1p19q_codeletion_status']=='Codel'].index
    non_codel = clin[clin['1p19q_codeletion_status']=='Non-codel'].index
    if len(codel)>=5 and len(non_codel)>=5:
        u, p_codel = stats.mannwhitneyu(zp3[codel], zp3[non_codel])
        print(f"1p/19q 共缺失 vs 非共缺失: p={p_codel:.4f}")
    
    # 生存分析
    print(f"\n  CGGA-325 生存分析:")
    surv = pd.DataFrame({'time': clin['OS'].values, 'event': clin['Censor (alive=0; dead=1)'].values, 'zp3': zp3.values}).dropna()
    median = surv['zp3'].median()
    high = surv[surv['zp3']>=median]
    low = surv[surv['zp3']<median]
    lr = logrank_test(high['time'], low['time'], event_observed_A=high['event'], event_observed_B=low['event'])
    print(f"  Log-rank(ZP3高/低): p={lr.p_value:.4f}")
    
    # Cox
    try:
        cph = CoxPHFitter()
        cph.fit(surv[['time','event','zp3']], duration_col='time', event_col='event')
        print(f"  Cox(连续): HR={cph.hazard_ratios_['zp3']:.3f}, p={cph.summary['p']['zp3']:.4f}")
    except Exception as e:
        print(f"  Cox失败: {e}")
    
    # 绘制对比图
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('CGGA-693 vs CGGA-325: ZP3-免疫特征相关性对比', fontsize=14, fontweight='bold')
    
    for ax, r, title in [(axes[0], r693, 'CGGA-693 (n=693)'), (axes[1], r325, 'CGGA-325 (n=325)')]:
        names = list(r['immune_scores'].keys())
        rhos = []
        ps = []
        for nm in names:
            rho, p = stats.spearmanr(r['zp3'], r['immune_scores'][nm][r['clinical'].index])
            rhos.append(rho)
            ps.append(p)
        colors = ['crimson' if p<0.05 else 'gray' for p in ps]
        ax.barh(names, rhos, color=colors, alpha=0.8, edgecolor='black')
        ax.axvline(0, color='black', linewidth=1)
        ax.set_xlim(-0.1, 0.6)
        ax.set_title(f'{title}')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fig_cgga693_vs_325_immune.png'), dpi=300, bbox_inches='tight')
    print(f"\n对比图保存: {os.path.join(output_dir, 'fig_cgga693_vs_325_immune.png')}")

if __name__ == '__main__':
    main()