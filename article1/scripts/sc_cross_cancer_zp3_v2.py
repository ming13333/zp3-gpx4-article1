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
Pan-cancer single-cell ZP3 validation v2 (Hypothesis 3: tissue specificity) — 5-cancer version
GBM + breast cancer + kidney cancer + melanoma + lung adenocarcinoma

v2 additions:
- Melanoma: Dissecting novel myeloid-derived cell states (pan-cancer myeloid atlas, 2.6GB)
  → Need to filter melanoma cells by disease column
- Lung adenocarcinoma: HTAN MSK LUAD (333MB) → Need to filter lung adenocarcinoma by disease
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
    """Find the first matching column in obs"""
    for c in candidates:
        if c in adata.obs.columns:
            return c
    return None

def get_gene_expr(adata, gene_key):
    """Extract gene expression (supports symbol or Ensembl)"""
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
    """Filter cells by disease keywords (for multi-disease mixed datasets)"""
    disease_col = find_col(adata, ['disease', 'disease_ontology_term_id', 'disease_ontology'])
    if disease_col is None:
        return adata, disease_col, None
    diseases = adata.obs[disease_col].astype(str)
    mask = pd.Series(False, index=adata.obs_names)
    for kw in keywords:
        mask = mask | diseases.str.lower().str.contains(kw, na=False)
    n_keep = mask.sum()
    print(f"    Disease column: {disease_col}, matching {kw} cells: {n_keep}/{adata.n_obs}")
    if n_keep == 0:
        return adata, disease_col, mask
    return adata, disease_col, mask

def analyze_dataset(adata, name, celltype_col, subset_mask=None):
    """Analyze a single dataset"""
    print(f"\n=== {name} ===")
    print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")
    
    # Apply disease filtering
    if subset_mask is not None:
        adata = adata[subset_mask]
        print(f"  Cells after filtering: {adata.n_obs}")
    
    zp3 = get_gene_expr(adata, 'ZP3')
    if zp3 is None:
        print("  ZP3 not found, skipping")
        return None
    
    print(f"  ZP3 detection rate (>0): {100*(zp3>0).mean():.2f}%")
    print(f"  ZP3 mean: {zp3.mean():.4f}, median: {np.median(zp3):.4f}")
    if (zp3 > 0).any():
        print(f"  ZP3 non-zero mean: {zp3[zp3>0].mean():.4f}")
    
    if celltype_col is None:
        print("  No cell type annotation, skipping")
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
    print(f"  ZP3+ cell type distribution (Top 5):")
    for _, row in ct_df.head(5).iterrows():
        print(f"    {str(row['cell_type'])[:40]}: {row['zp3_pos_pct']:.2f}% (n={row['n_cells']})")
    
    # Immune marker co-expression
    print(f"  Immune marker co-expression (ZP3+ vs ZP3-):")
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
                    print(f"    {gene}: ZP3+ {pc_pos:.1f}% vs ZP3- {pc_neg:.1f}% (OR={oddsratio:.2f}, p={p:.2e})")
                    coexpr_rows.append({'dataset': name, 'gene': gene, 'or': oddsratio, 'p': p,
                                        'zp3pos_pct': pc_pos, 'zp3neg_pct': pc_neg})
                except Exception:
                    pass
    
    return ct_df, coexpr_rows

def main():
    print("=" * 60)
    print("Cross-cancer single-cell ZP3 validation v2 (5 cancer types)")
    print("=" * 60)
    
    all_ct = []
    all_coexpr = []
    
    # 1. GBM
    print("\n[1] GBM (GSE141982)...")
    try:
        gbm = sc.read_h5ad(os.path.join(ROOT, "output", "h1_pilot", "h1_adata.h5ad"))
        ct_col = find_col(gbm, ['cell_type', 'celltype', 'labels_unif'])
        print(f"  Cell type column: {ct_col}")
        r = analyze_dataset(gbm, "GBM (GSE141982)", ct_col)
        if r: all_ct.append(r[0]); all_coexpr.extend(r[1])
    except Exception as e:
        print(f"  Failed: {e}")
    
    # 2. Breast Cancer
    print("\n[2] Breast Cancer (HTAN)...")
    bp = os.path.join(SC_DATA_DIR, "breast_htan.h5ad")
    if os.path.exists(bp):
        try:
            b = sc.read_h5ad(bp)
            ct_col = find_col(b, ['cell_type', 'celltype', 'labels_unif'])
            print(f"  Cell type column: {ct_col}")
            r = analyze_dataset(b, "Breast Cancer (HTAN)", ct_col)
            if r: all_ct.append(r[0]); all_coexpr.extend(r[1])
        except Exception as e:
            print(f"  Failed: {e}")
    else:
        print("  Data not found")
    
    # 3. Renal Cell Carcinoma
    print("\n[3] Renal Cell Carcinoma (HTAN)...")
    rp = os.path.join(SC_DATA_DIR, "rcc_htan.h5ad")
    if os.path.exists(rp):
        try:
            rcc = sc.read_h5ad(rp)
            ct_col = find_col(rcc, ['cell_type', 'celltype', 'labels_unif'])
            print(f"  Cell type column: {ct_col}")
            r = analyze_dataset(rcc, "Renal Cell Carcinoma (HTAN)", ct_col)
            if r: all_ct.append(r[0]); all_coexpr.extend(r[1])
        except Exception as e:
            print(f"  Failed: {e}")
    else:
        print("  Data not found")
    
    # 4. Melanoma (multi-cancer dataset, filter melanoma)
    print("\n[4] Melanoma (myeloid atlas)...")
    mp = os.path.join(SC_DATA_DIR, "melanoma_myeloid.h5ad")
    if os.path.exists(mp):
        try:
            mel = sc.read_h5ad(mp)
            print(f"  Total cells: {mel.n_obs}, obs columns: {list(mel.obs.columns)[:12]}")
            ct_col = find_col(mel, ['cell_type', 'celltype', 'labels_unif', 'compartment'])
            print(f"  Cell type column: {ct_col}")
            mel2, dcol, mask = filter_by_disease(mel, ['melanoma'])
            if dcol and mask is not None and mask.sum() > 0:
                r = analyze_dataset(mel2, "Melanoma", ct_col, mask)
            else:
                r = analyze_dataset(mel2, "Melanoma", ct_col)
            if r: all_ct.append(r[0]); all_coexpr.extend(r[1])
        except Exception as e:
            print(f"  Failed: {e}")
    else:
        print("  Data not found (downloading)")
    
    # 5. Lung adenocarcinoma
    print("\n[5] Lung Adenocarcinoma (HTAN MSK)...")
    lp = os.path.join(SC_DATA_DIR, "luad_htan.h5ad")
    if os.path.exists(lp):
        try:
            luad = sc.read_h5ad(lp)
            print(f"  Total cells: {luad.n_obs}, obs columns: {list(luad.obs.columns)[:12]}")
            ct_col = find_col(luad, ['cell_type', 'celltype', 'labels_unif'])
            print(f"  Cell type column: {ct_col}")
            luad2, dcol, mask = filter_by_disease(luad, ['lung adenocarcinoma'])
            if dcol and mask is not None and mask.sum() > 0:
                r = analyze_dataset(luad2, "Lung Adenocarcinoma (HTAN MSK)", ct_col, mask)
            else:
                r = analyze_dataset(luad2, "Lung Adenocarcinoma (HTAN MSK)", ct_col)
            if r: all_ct.append(r[0]); all_coexpr.extend(r[1])
        except Exception as e:
            print(f"  Failed: {e}")
    else:
        print("  Data not found (downloading)")
    
    # Summary
    if all_ct:
        combined = pd.concat(all_ct, ignore_index=True)
        combined.to_csv(os.path.join(OUTPUT_DIR, "sc_cross_cancer_zp3_celltype_v2.csv"), index=False)
        print(f"\nSaved: sc_cross_cancer_zp3_celltype_v2.csv")
        
        if all_coexpr:
            coexpr_df = pd.DataFrame(all_coexpr)
            coexpr_df.to_csv(os.path.join(OUTPUT_DIR, "sc_cross_cancer_zp3_coexpr_v2.csv"), index=False)
            print(f"Saved: sc_cross_cancer_zp3_coexpr_v2.csv")
        
        visualize(combined, all_coexpr)
    else:
        print("\nNo data available!")
    
    print("\n" + "=" * 60)
    print("Cross-cancer single-cell validation v2 complete")
    print("=" * 60)

def visualize(combined, all_coexpr):
    """Visualize cross-cancer results"""
    plt.style.use('seaborn-v0_8-whitegrid')
    
    n_panels = 2 + (1 if all_coexpr else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(6*n_panels, 5.5))
    
    # 1. Cell type with highest ZP3+ proportion per cancer
    ax = axes[0]
    ds_max = combined.groupby('dataset').apply(lambda x: x.loc[x['zp3_pos_pct'].idxmax()])
    ax.barh(ds_max['dataset'], ds_max['zp3_pos_pct'], color='#d62728')
    for i, (_, row) in enumerate(ds_max.iterrows()):
        ax.text(row['zp3_pos_pct'] + 0.1, i, f"{str(row['cell_type'])[:18]} {row['zp3_pos_pct']:.1f}%", 
                va='center', fontsize=9)
    ax.set_xlabel('ZP3+ proportion (%)')
    ax.set_title('Top ZP3+ cell type by cancer')
    
    # 2. Myeloid vs non-myeloid ZP3+ proportion
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
    
    # 3. Coexpression OR
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
    print(f"Saved: {fig_path}")

if __name__ == "__main__":
    main()
