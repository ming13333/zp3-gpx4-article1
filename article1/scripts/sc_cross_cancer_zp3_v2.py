#!/usr/bin/env python
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
跨癌种单细胞 ZP3 验证 v2（假说 3：组织特异性）— 5 癌种版
GBM + 乳腺癌 + 肾癌 + 黑色素瘤 + 肺腺癌

v2 新增：
- 黑色素瘤：Dissecting novel myeloid-derived cell states (多癌种髓系图谱, 2.6GB)
  → 需按 disease 列筛选 melanoma 细胞
- 肺腺癌：HTAN MSK LUAD (333MB) → 需按 disease 筛选 lung adenocarcinoma
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import sys
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR = os.path.join(ROOT, "article1", "results", "sc_cross_cancer")
SC_DATA_DIR = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling", "sc_data")
ZP3_ENSEMBL = "ENSG00000188372"

IMMUNE_MARKERS = ['TREM2', 'CD68', 'CD163', 'MRC1', 'CSF1R']

def find_col(adata, candidates):
    """在 obs 中找第一个匹配的列"""
    for c in candidates:
        if c in adata.obs.columns:
            return c
    return None

def get_gene_expr(adata, gene_key):
    """提取基因表达（支持 symbol 或 Ensembl）"""
    if gene_key in adata.var_names:
        idx = list(adata.var_names).index(gene_key)
    elif ZP3_ENSEMBL in adata.var_names and gene_key == 'ZP3':
        idx = list(adata.var_names).index(ZP3_ENSEMBL)
    else:
        return None
    X = adata.X[:, idx]
    if hasattr(X, 'toarray'):
        X = X.toarray().flatten()
    else:
        X = np.asarray(X).flatten()
    return X

def filter_by_disease(adata, keywords):
    """按疾病关键词筛选细胞（用于多疾病混合数据集）"""
    disease_col = find_col(adata, ['disease', 'disease_ontology_term_id', 'disease_ontology'])
    if disease_col is None:
        return adata, disease_col, None
    diseases = adata.obs[disease_col].astype(str)
    mask = pd.Series(False, index=adata.obs_names)
    for kw in keywords:
        mask = mask | diseases.str.lower().str.contains(kw, na=False)
    n_keep = mask.sum()
    print(f"    疾病列: {disease_col}, 匹配 {kw} 细胞: {n_keep}/{adata.n_obs}")
    if n_keep == 0:
        return adata, disease_col, mask
    return adata, disease_col, mask

def analyze_dataset(adata, name, celltype_col, subset_mask=None):
    """分析单个数据集"""
    print(f"\n=== {name} ===")
    print(f"  细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    
    # 应用疾病筛选
    if subset_mask is not None:
        adata = adata[subset_mask]
        print(f"  筛选后细胞数: {adata.n_obs}")
    
    zp3 = get_gene_expr(adata, 'ZP3')
    if zp3 is None:
        print("  ZP3 未找到，跳过")
        return None
    
    print(f"  ZP3 检出率 (>0): {100*(zp3>0).mean():.2f}%")
    print(f"  ZP3 均值: {zp3.mean():.4f}, 中位数: {np.median(zp3):.4f}")
    if (zp3 > 0).any():
        print(f"  ZP3 非零值均值: {zp3[zp3>0].mean():.4f}")
    
    if celltype_col is None:
        print("  无细胞类型注释，跳过")
        return None
    
    results = []
    celltypes = adata.obs[celltype_col].astype(str)
    
    for ct in celltypes.unique():
        mask = celltypes == ct
        n_cells = mask.sum()
        ct_zp3 = zp3[mask]
        n_pos = (ct_zp3 > 0).sum()
        results.append({
            'dataset': name,
            'cell_type': ct,
            'n_cells': n_cells,
            'zp3_pos_pct': 100 * n_pos / n_cells if n_cells > 0 else 0,
            'zp3_mean': ct_zp3.mean() if n_cells > 0 else 0,
            'zp3_pos_mean': ct_zp3[ct_zp3 > 0].mean() if n_pos > 0 else 0
        })
    
    ct_df = pd.DataFrame(results).sort_values('zp3_pos_pct', ascending=False)
    print(f"  ZP3+ 细胞类型分布 (Top 5):")
    for _, row in ct_df.head(5).iterrows():
        print(f"    {str(row['cell_type'])[:40]}: {row['zp3_pos_pct']:.2f}% (n={row['n_cells']})")
    
    # 免疫标志物共表达
    print(f"  免疫标志物共表达 (ZP3+ vs ZP3-):")
    coexpr_rows = []
    for gene in IMMUNE_MARKERS:
        gene_expr = get_gene_expr(adata, gene)
        if gene_expr is None:
            continue
        zp3_pos = zp3 > 0
        gene_pos = gene_expr > 0
        if zp3_pos.sum() > 0 and gene_pos.sum() > 0:
            a = (zp3_pos & gene_pos).sum()
            b = (zp3_pos & ~gene_pos).sum()
            c = (~zp3_pos & gene_pos).sum()
            d = (~zp3_pos & ~gene_pos).sum()
            if a + b > 0 and c + d > 0:
                try:
                    oddsratio, p = stats.fisher_exact([[a, b], [c, d]])
                    pc_pos = 100 * a / (a + b)
                    pc_neg = 100 * c / (c + d)
                    print(f"    {gene}: ZP3+中{pc_pos:.1f}% vs ZP3-中{pc_neg:.1f}% (OR={oddsratio:.2f}, p={p:.2e})")
                    coexpr_rows.append({'dataset': name, 'gene': gene, 'or': oddsratio, 'p': p,
                                        'zp3pos_pct': pc_pos, 'zp3neg_pct': pc_neg})
                except Exception:
                    pass
    
    return ct_df, coexpr_rows

def main():
    print("=" * 60)
    print("跨癌种单细胞 ZP3 验证 v2 (5 癌种)")
    print("=" * 60)
    
    all_ct = []
    all_coexpr = []
    
    # 1. GBM
    print("\n[1] GBM (GSE141982)...")
    try:
        gbm = sc.read_h5ad(os.path.join(ROOT, "output", "h1_pilot", "h1_adata.h5ad"))
        ct_col = find_col(gbm, ['cell_type', 'celltype', 'labels_unif'])
        print(f"  细胞类型列: {ct_col}")
        r = analyze_dataset(gbm, "GBM (GSE141982)", ct_col)
        if r: all_ct.append(r[0]); all_coexpr.extend(r[1])
    except Exception as e:
        print(f"  失败: {e}")
    
    # 2. 乳腺癌
    print("\n[2] Breast Cancer (HTAN)...")
    bp = os.path.join(SC_DATA_DIR, "breast_htan.h5ad")
    if os.path.exists(bp):
        try:
            b = sc.read_h5ad(bp)
            ct_col = find_col(b, ['cell_type', 'celltype', 'labels_unif'])
            print(f"  细胞类型列: {ct_col}")
            r = analyze_dataset(b, "Breast Cancer (HTAN)", ct_col)
            if r: all_ct.append(r[0]); all_coexpr.extend(r[1])
        except Exception as e:
            print(f"  失败: {e}")
    else:
        print("  数据不存在")
    
    # 3. 肾癌
    print("\n[3] Renal Cell Carcinoma (HTAN)...")
    rp = os.path.join(SC_DATA_DIR, "rcc_htan.h5ad")
    if os.path.exists(rp):
        try:
            rcc = sc.read_h5ad(rp)
            ct_col = find_col(rcc, ['cell_type', 'celltype', 'labels_unif'])
            print(f"  细胞类型列: {ct_col}")
            r = analyze_dataset(rcc, "Renal Cell Carcinoma (HTAN)", ct_col)
            if r: all_ct.append(r[0]); all_coexpr.extend(r[1])
        except Exception as e:
            print(f"  失败: {e}")
    else:
        print("  数据不存在")
    
    # 4. 黑色素瘤（多癌种数据集，筛选 melanoma）
    print("\n[4] Melanoma (myeloid atlas)...")
    mp = os.path.join(SC_DATA_DIR, "melanoma_myeloid.h5ad")
    if os.path.exists(mp):
        try:
            mel = sc.read_h5ad(mp)
            print(f"  总细胞数: {mel.n_obs}, obs 列: {list(mel.obs.columns)[:12]}")
            ct_col = find_col(mel, ['cell_type', 'celltype', 'labels_unif', 'compartment'])
            print(f"  细胞类型列: {ct_col}")
            mel2, dcol, mask = filter_by_disease(mel, ['melanoma'])
            if dcol and mask is not None and mask.sum() > 0:
                r = analyze_dataset(mel2, "Melanoma", ct_col, mask)
            else:
                r = analyze_dataset(mel2, "Melanoma", ct_col)
            if r: all_ct.append(r[0]); all_coexpr.extend(r[1])
        except Exception as e:
            print(f"  失败: {e}")
    else:
        print("  数据不存在（下载中）")
    
    # 5. 肺腺癌
    print("\n[5] Lung Adenocarcinoma (HTAN MSK)...")
    lp = os.path.join(SC_DATA_DIR, "luad_htan.h5ad")
    if os.path.exists(lp):
        try:
            luad = sc.read_h5ad(lp)
            print(f"  总细胞数: {luad.n_obs}, obs 列: {list(luad.obs.columns)[:12]}")
            ct_col = find_col(luad, ['cell_type', 'celltype', 'labels_unif'])
            print(f"  细胞类型列: {ct_col}")
            luad2, dcol, mask = filter_by_disease(luad, ['lung adenocarcinoma'])
            if dcol and mask is not None and mask.sum() > 0:
                r = analyze_dataset(luad2, "Lung Adenocarcinoma (HTAN MSK)", ct_col, mask)
            else:
                r = analyze_dataset(luad2, "Lung Adenocarcinoma (HTAN MSK)", ct_col)
            if r: all_ct.append(r[0]); all_coexpr.extend(r[1])
        except Exception as e:
            print(f"  失败: {e}")
    else:
        print("  数据不存在（下载中）")
    
    # 汇总
    if all_ct:
        combined = pd.concat(all_ct, ignore_index=True)
        combined.to_csv(os.path.join(OUTPUT_DIR, "sc_cross_cancer_zp3_celltype_v2.csv"), index=False)
        print(f"\n已保存: sc_cross_cancer_zp3_celltype_v2.csv")
        
        if all_coexpr:
            coexpr_df = pd.DataFrame(all_coexpr)
            coexpr_df.to_csv(os.path.join(OUTPUT_DIR, "sc_cross_cancer_zp3_coexpr_v2.csv"), index=False)
            print(f"已保存: sc_cross_cancer_zp3_coexpr_v2.csv")
        
        visualize(combined, all_coexpr)
    else:
        print("\n无可用数据！")
    
    print("\n" + "=" * 60)
    print("跨癌种单细胞验证 v2 完成")
    print("=" * 60)

def visualize(combined, all_coexpr):
    """可视化跨癌种结果"""
    plt.style.use('seaborn-v0_8-whitegrid')
    
    n_panels = 2 + (1 if all_coexpr else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(6*n_panels, 5.5))
    
    # 1. 每癌种 ZP3+ 比例最高的细胞类型
    ax = axes[0]
    ds_max = combined.groupby('dataset').apply(lambda x: x.loc[x['zp3_pos_pct'].idxmax()])
    ax.barh(ds_max['dataset'], ds_max['zp3_pos_pct'], color='#d62728')
    for i, (_, row) in enumerate(ds_max.iterrows()):
        ax.text(row['zp3_pos_pct'] + 0.1, i, f"{str(row['cell_type'])[:18]} {row['zp3_pos_pct']:.1f}%", 
                va='center', fontsize=9)
    ax.set_xlabel('ZP3+ proportion (%)')
    ax.set_title('Top ZP3+ cell type by cancer')
    
    # 2. 髓系 vs 非髓系 ZP3+ 比例
    ax = axes[1]
    mye_kw = ['macro', 'mono', 'myeloid', 'tam', 'mg', 'microglia', 'dendritic', 'dc ']
    rows = []
    for _, row in combined.iterrows():
        is_mye = any(k in str(row['cell_type']).lower() for k in mye_kw)
        rows.append({'dataset': row['dataset'], 'compartment': 'Myeloid' if is_mye else 'Other', 
                     'zp3_pos_pct': row['zp3_pos_pct']})
    mdf = pd.DataFrame(rows)
    pivot = mdf.groupby(['dataset', 'compartment'])['zp3_pos_pct'].max().unstack(fill_value=0)
    pivot = pivot.reindex(sorted(pivot.index))
    pivot.plot(kind='bar', ax=ax, color=['#1f77b4', '#ff7f0e'])
    ax.set_ylabel('Max ZP3+ proportion (%)')
    ax.set_title('Myeloid vs Other cell types')
    ax.legend(title='')
    ax.tick_params(axis='x', rotation=30)
    
    # 3. 共表达 OR
    if all_coexpr:
        ax = axes[2]
        coexpr_df = pd.DataFrame(all_coexpr)
        plot_df = coexpr_df[coexpr_df['gene'] == 'TREM2'].copy()
        if len(plot_df) == 0:
            plot_df = coexpr_df[coexpr_df['gene'].isin(['TREM2', 'CSF1R'])].copy()
        if len(plot_df) > 0:
            plot_df = plot_df.sort_values('or')
            ax.barh(plot_df['dataset'], plot_df['or'], color='#9467bd')
            ax.set_xlabel('Fisher OR (log scale)')
            ax.set_xscale('log')
            ax.set_title('ZP3+ vs ZP3- TREM2 coexpression OR')
        else:
            ax.text(0.5, 0.5, 'No TREM2 data', ha='center', va='center')
            ax.axis('off')
    
    plt.tight_layout()
    fig_path = os.path.join(OUTPUT_DIR, "fig_sc_cross_cancer_zp3_v2.png")
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"已保存: {fig_path}")

if __name__ == "__main__":
    main()
