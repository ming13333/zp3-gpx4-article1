#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
IMvigor210 Cohort ZP3 and Immunotherapy Response Analysis
Analyzing the predictive value of ZP3 in urothelial cancer immunotherapy
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import mannwhitneyu, chi2_contingency
import warnings
warnings.filterwarnings('ignore')

# Set Chinese font display
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

print("=" * 60)
print("IMvigor210 Cohort ZP3 Immunotherapy Response Analysis")
print("=" * 60)

# ============================================================
# 1. Load data
# ============================================================
print("\n[1] Loading data...")

# Load phenotype data
IMDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'output', 'immunotherapy_validation')
pData = pd.read_csv(os.path.join(IMDIR, 'pData_IMvigor210.csv'))
print(f"  Phenotype data: {pData.shape[0]} samples, {pData.shape[1]} variables")

# Load feature data (gene annotations)
fData = pd.read_csv(os.path.join(IMDIR, 'fData_IMvigor210.csv'))
print(f"  Feature data: {fData.shape[0]} genes")

# Load expression matrix
print("  Loading expression matrix (this may take some time)...")
exmat = pd.read_csv(os.path.join(IMDIR, 'exmat_censored_IMvigor210.csv'), index_col=0)
print(f"  Expression matrix: {exmat.shape[0]} genes × {exmat.shape[1]} samples")

# ============ Explicit row-matching validation (2026-08-17 scientific audit fix) ============
# exmat row names are gene_1..N (written in fData order). Must verify before analysis:
#  (a) row count matches fData; (b) the ZP3 row (fData Symbol column) can be located and is reasonably expressed.
assert exmat.shape[0] == fData.shape[0], \
    f"Row count mismatch: exmat {exmat.shape[0]} vs fData {fData.shape[0]}"
fData_tmp = fData.copy()
fData_tmp['row_name'] = ['gene_' + str(i) for i in range(1, len(fData_tmp) + 1)]
_sym2row = dict(zip(fData_tmp['Symbol'].astype(str), fData_tmp['row_name']))
_zp3_row = _sym2row.get('ZP3')
assert _zp3_row is not None, "ZP3 symbol not found in fData"
_zp3_val = exmat.loc[_zp3_row].astype(float)
assert (_zp3_val > 0).mean() > 0.5, f"ZP3 row non-zero proportion too low: {(_zp3_val > 0).mean():.1%}"
print(f"  [Row-matching validation passed] fData rows={fData.shape[0]}, ZP3 located at row '{_zp3_row}' (Entrez 7784), "
      f"non-zero proportion={(_zp3_val > 0).mean():.1%}, median={_zp3_val.median():.1f}")

# ============================================================
# 2. Data preprocessing
# ============================================================
print("\n[2] Data preprocessing...")

# Create gene name mapping (from fData)
gene_map = fData.set_index('entrez_id')['symbol'].to_dict()

# Convert expression matrix row names from gene_XXX to actual gene names
# Expression matrix row names are gene_1, gene_2... corresponding to the order of fData
exmat.index = [gene_map.get(fData.iloc[i, 0], f"gene_{i}") for i in range(len(exmat))]
print(f"  Expression matrix row names updated to gene symbols")

# Extract ZP3 expression
if 'ZP3' in exmat.index:
    zp3_expr = exmat.loc['ZP3']
    print(f"  ZP3 expression extracted: {len(zp3_expr)} samples")
    print(f"  ZP3 expression range: {zp3_expr.min():.2f} - {zp3_expr.max():.2f}")
    print(f"  ZP3 median expression: {zp3_expr.median():.2f}")
else:
    print("  Error: ZP3 gene not found!")
    exit(1)

# Merge clinical information
pData_aligned = pData.set_index('X').loc[exmat.columns]
print(f"  Number of samples after alignment: {len(pData_aligned)}")

# ============================================================
# 3. ZP3 expression and treatment response analysis
# ============================================================
print("\n[3] ZP3 expression and treatment response analysis...")

# Add ZP3 expression to clinical data
pData_aligned['ZP3_expr'] = zp3_expr.values

# Define response groups (exclude NE - not evaluable)
responders = pData_aligned[pData_aligned['binaryResponse'] == 'CR/PR']
non_responders = pData_aligned[pData_aligned['binaryResponse'] == 'SD/PD']

print(f"  Responders (CR/PR): {len(responders)} cases")
print(f"  Non-responders (SD/PD): {len(non_responders)} cases")

# Difference in ZP3 expression between response groups
zp3_resp = responders['ZP3_expr']
zp3_nonresp = non_responders['ZP3_expr']

# Mann-Whitney U test
stat_u, p_u = mannwhitneyu(zp3_resp, zp3_nonresp, alternative='two-sided')
print(f"\n  ZP3 expression difference test (Mann-Whitney U):")
print(f"    Responders median: {zp3_resp.median():.2f} (IQR: {zp3_resp.quantile(0.25):.2f} - {zp3_resp.quantile(0.75):.2f})")
print(f"    Non-responders median: {zp3_nonresp.median():.2f} (IQR: {zp3_nonresp.quantile(0.25):.2f} - {zp3_nonresp.quantile(0.75):.2f})")
print(f"    U statistic: {stat_u:.2f}, p-value: {p_u:.4f}")

# Effect size (rank-biserial correlation)
# 2026-08-17 audit fix: unify r = 2U/(n1*n2) - 1 (U for responder group), positive = responders higher
n1, n2 = len(zp3_resp), len(zp3_nonresp)
r = (2 * stat_u) / (n1 * n2) - 1
print(f"    Effect size (rank-biserial r): {r:.4f}  (positive = higher ZP3 in responders)")

# ============================================================
# 4. ZP3 high/low expression group analysis
# ============================================================
print("\n[4] ZP3 high/low expression group analysis...")

# Split using median
median_zp3 = pData_aligned['ZP3_expr'].median()
pData_aligned['ZP3_group'] = pData_aligned['ZP3_expr'].apply(lambda x: 'High' if x >= median_zp3 else 'Low')

# Contingency table
contingency = pd.crosstab(pData_aligned['ZP3_group'], pData_aligned['binaryResponse'])
print("\n  ZP3 group × treatment response contingency table:")
print(contingency)

# Chi-square test
chi2, p_chi, dof, expected = chi2_contingency(contingency)
print(f"\n  Chi-square test: χ² = {chi2:.2f}, p = {p_chi:.4f}")

# Calculate response rates
response_rates = pData_aligned.groupby('ZP3_group')['binaryResponse'].apply(
    lambda x: (x == 'CR/PR').sum() / len(x) * 100
)
print(f"\n  Response rates by group:")
for group, rate in response_rates.items():
    print(f"    {group} ZP3: {rate:.1f}%")

# ============================================================
# 5. Survival analysis
# ============================================================
print("\n[5] Survival analysis...")

# Check survival data
if 'OS' in pData_aligned.columns and 'censored' in pData_aligned.columns:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    # KM curves
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

    # Log-rank test
    results = logrank_test(pData_aligned.loc[mask_high, 'OS'],
                           pData_aligned.loc[mask_low, 'OS'],
                           pData_aligned.loc[mask_high, 'censored'],
                           pData_aligned.loc[mask_low, 'censored'])

    print(f"  Median survival:")
    print(f"    ZP3 High: {kmf_high.median_survival_time_:.1f} days")
    print(f"    ZP3 Low: {kmf_low.median_survival_time_:.1f} days")
    print(f"  Log-rank test: χ² = {results.test_statistic:.2f}, p = {results.p_value:.4f}")

    # Plot KM curves
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
    print("  KM curve saved: fig_imvigor210_km_zp3.png")
else:
    print("  Warning: Survival data columns incomplete, skipping survival analysis")

# ============================================================
# 6. Immune signature correlation analysis
# ============================================================
print("\n[6] Immune signature correlation analysis...")

# Define immune signature gene sets
immune_signatures = {
    'T_cell_exhaustion': ['PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'TOX'],
    'Cytolytic_activity': ['GZMA', 'GZMB', 'PRF1', 'IFNG'],
    'Immune_checkpoint': ['CD274', 'PDCD1LG2', 'CTLA4', 'PDCD1', 'LAG3', 'HAVCR2'],
    'Macrophage_M2': ['CD163', 'MSR1', 'MRC1', 'VSIG4', 'CD206'],
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF18'],
    'TGF_beta': ['TGFB1', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2'],
    'IFN_gamma_response': ['STAT1', 'IRF1', 'CXCL10', 'CXCL9', 'IDO1']
}

# Calculate mean expression for each signature
print("\n  Correlation between immune signatures and ZP3:")
correlation_results = []

for sig_name, genes in immune_signatures.items():
    # Check whether genes exist
    valid_genes = [g for g in genes if g in exmat.index]
    if len(valid_genes) >= 2:
        # Calculate signature score (mean expression)
        sig_score = exmat.loc[valid_genes].mean(axis=0)
        # Compute Spearman correlation with ZP3
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
# 7. Visualization
# ============================================================
print("\n[7] Generating visualizations...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 7.1 ZP3 expression distribution (by response group)
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

# 7.2 Response rate in ZP3 high/low groups
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

# 7.3 Immune signature correlation heatmap
ax3 = axes[1, 0]
if correlation_results:
    corr_df = pd.DataFrame(correlation_results)
    corr_df = corr_df.set_index('Signature')
    sns.heatmap(corr_df[['Spearman_rho']], annot=True, fmt='.3f', cmap='RdBu_r',
                center=0, ax=ax3, cbar_kws={'label': 'Spearman ρ'})
    ax3.set_title('ZP3 vs Immune Signatures', fontsize=12)
    ax3.set_ylabel('')

# 7.4 ZP3 expression distribution (overall)
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
print("  Main figure saved: fig_imvigor210_zp3_analysis.png")

# ============================================================
# 8. Save results
# ============================================================
print("\n[8] Saving results...")

# Save detailed results
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

# Save as CSV
results_df = pd.DataFrame([results_summary])
results_df.to_csv(os.path.join(IMDIR, 'imvigor210_zp3_results.csv'), index=False)
print("  Results saved: imvigor210_zp3_results.csv")

# Save immune signature correlations
if correlation_results:
    corr_df = pd.DataFrame(correlation_results)
    corr_df.to_csv(os.path.join(IMDIR, 'imvigor210_zp3_immune_correlations.csv'), index=False)
    print("  Immune correlations saved: imvigor210_zp3_immune_correlations.csv")

# ============================================================
# 9. Summary
# ============================================================
print("\n" + "=" * 60)
print("Analysis Summary")
print("=" * 60)

print(f"""
Cohort: IMvigor210 (Urothelial carcinoma, Atezolizumab treatment)
Sample size: {len(pData_aligned)} cases (CR/PR: {len(responders)}, SD/PD: {len(non_responders)})

Key findings:
1. ZP3 expression difference:
   - Responders vs non-responders: p = {p_u:.4f}
   - Effect size r = {r:.4f}

2. Treatment response rate:
   - ZP3 Low: {response_rates.get('Low', 0):.1f}%
   - ZP3 High: {response_rates.get('High', 0):.1f}%
   - Chi-square test p = {p_chi:.4f}

Conclusion:
{"ZP3 high expression is associated with poorer immunotherapy response" if r > 0 and p_u < 0.05 else
 "ZP3 low expression is associated with poorer immunotherapy response" if r < 0 and p_u < 0.05 else
 "ZP3 expression is not significantly associated with immunotherapy response"}
""")

print("\nAnalysis complete!")
