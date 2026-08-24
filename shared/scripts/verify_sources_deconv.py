# -*- coding: utf-8 -*-
"""
Method 1 (single-cell deconvolution/source identification) + response to the reviewer's concern that "bulk ZP3 source = tumor vs myeloid contamination".

Core idea (not relying on fragile cell-type reference profile deconvolution, instead using two layers of targeted tests to directly address the confounding concern):

  [Test A] Single-cell level ZP3↔TREM2 co-enrichment (local h5ad, zero network)
     - Direct test: whether ZP3+ cells are enriched for TREM2+ (all cells / within myeloid / within MG-TAM-DC subsets)
     - If ZP3+ cells within myeloid almost all co-express TREM2, then the "bulk ZP3↔TREM2 correlation" is not
       caused by tumor cell contamination, but rather the same population of TREM2+ myeloid cells simultaneously contributes ZP3 and TREM2.
     - Use Fisher's exact test + enrichment strength (OR) + univariate logistic regression (logistic OR per unit ZP3)

  [Test B] bulk level: whether the ZP3-TREM2 association is independent of “total myeloid burden”
     - Use pan-myeloid markers (CD68/CD14/LYZ/CSF1R/ITGAM) to construct a myeloid index
     - Pearson: ZP3 vs myeloid index
     - Partial correlation: ZP3 vs TREM2 after controlling for the myeloid index (implemented via residual correlation)
     - If the association attenuates substantially after controlling for myeloid burden -> confounding (the apparent bulk-level ZP3↔TREM2 correlation comes from total myeloid abundance)
     - If attenuated but still significant -> there is a specific association independent of total myeloid abundance
     - Both results are honest \"publishable\" conclusions, but the interpretation differs.

  2026-08-10  Craft mode execution
"""
import numpy as np
import pandas as pd
import scipy.stats as stats
import anndata
import os

BASE = os.path.dirname(os.path.abspath(__file__))

def _project_root():
    d = os.path.dirname(os.path.abspath(__file__))
    while True:
        if os.path.isdir(os.path.join(d, "output")):
            return d
        p = os.path.dirname(d)
        if p == d:
            break
        d = p
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = _project_root()
OUT = []

def log(msg):
    print(msg)
    OUT.append(str(msg))

# ---------------------------------------------------------------
# [Test A] Single-cell-level ZP3 ↔ TREM2 co-enrichment
# ---------------------------------------------------------------
def main_A():
    h5ad_path = os.path.join(ROOT, "output", "h1_pilot", "h1_adata_subtyped.h5ad")
    adata = anndata.read_h5ad(h5ad_path)
    X = adata.X
    import scipy.sparse as sp
    if sp.issparse(X):
        zp3 = X[:, adata.var_names.get_loc('ZP3')].toarray().ravel()
        trem2 = X[:, adata.var_names.get_loc('TREM2')].toarray().ravel()
    else:
        zp3 = X[:, adata.var_names.get_loc('ZP3')].ravel()
        trem2 = X[:, adata.var_names.get_loc('TREM2')].ravel()
    obs = adata.obs
    zp3_pos = zp3 > 0
    trem2_pos = trem2 > 0
    log("=" * 70)
    log("[Test A] Single-cell level: whether ZP3+ cells are enriched for TREM2+  (n_cells=%d)" % len(zp3))
    log("  var space = %d HVG (non-whole-genome); ZP3/TREM2 both retained due to high variability" % adata.n_vars)
    log("  Background: ZP3+ %d cells (%.2f%%), TREM2+ %d cells (%.2f%%)"
        % (zp3_pos.sum(), 100*zp3_pos.mean(), trem2_pos.sum(), 100*trem2_pos.mean()))

    def fisher_or_table(a, b, label):
        a = np.asarray(a, bool); b = np.asarray(b, bool)
        aa = int(((a) & (b)).sum()); ab = int(((a) & (~b)).sum())
        ba = int(((~a) & (b)).sum()); bb = int(((~a) & (~b)).sum())
        OR, p = stats.fisher_exact([[aa, ab], [ba, bb]])
        frac_zp3_trem2 = aa/(aa+ab) if (aa+ab)>0 else float('nan')
        frac_bg_trem2 = (aa+ba)/len(a) if len(a)>0 else float('nan')
        log("  [%s]  ZP3+&TREM2+=%d, ZP3+&TREM2-=%d, ZP3-&TREM2+=%d, ZP3-&TREM2-=%d"
            % (label, aa, ab, ba, bb))
        log("      TREM2+ proportion in ZP3+=%.1f%% | overall background TREM2+ proportion=%.1f%%" % (100*frac_zp3_trem2, 100*frac_bg_trem2))
        log("      Fisher OR=%.2f, p=%.3g" % (OR, p))
        return OR, p

    # A0 all cells
    fisher_or_table(zp3_pos, trem2_pos, "all cells")
    # A1 within myeloid
    myeloid = obs['myeloid'].values.astype(bool)
    if myeloid.sum() > 0:
        m = fisher_or_table(zp3_pos[myeloid], trem2_pos[myeloid], "within myeloid myel=True")
    # A2 stratified by myeloid subpopulation
    sub = obs['myeloid_subclass'].values
    sub = np.array([str(x) for x in sub])
    for s in ['MG', 'TAM', 'DC']:
        sel = sub == s
        if sel.sum() >= 10:
            fp3 = zp3_pos[sel]; ftrem = trem2_pos[sel]
            aa=int((fp3&ftrem).sum()); ab=int((fp3&~ftrem).sum())
            ba=int((~fp3&ftrem).sum()); bb=int((~fp3&~ftrem).sum())
            if (aa+ab)>0 and (ba+bb)>0:
                OR,p = stats.fisher_exact([[aa,ab],[ba,bb]])
                log("  [subpopulation:%s n=%d]  ZP3+&TREM2+=%d/%d  | within-subpopulation TREM2+ background=%.1f%% | Fisher OR=%.2f p=%.3g"
                    % (s, sel.sum(), aa, aa+ab, 100*ftrem.mean(), OR, p))
    # A3 univariate logistic regression (ZP3 expression -> TREM2+, all cells) gives OR per unit ZP3
    try:
        from scipy.optimize import minimize
        y = trem2_pos.astype(float); x = zp3
        Xd = np.column_stack([np.ones_like(x), x])
        def negll(b):
            z = Xd @ b
            p = 1/(1+np.exp(-z))
            p = np.clip(p, 1e-12, 1-1e-12)
            return -(y*np.log(p) + (1-y)*np.log(1-p)).sum()
        res = minimize(negll, np.zeros(2), method='BFGS')
        b0, b1 = res.x
        log("  [logistic regression, unadjusted] logit(TREM2+)=%.3f+%.3f*ZP3 -> OR per +1(log1p) of ZP3=%.2f"
            % (b0, b1, np.exp(b1)))
    except Exception as e:
        log("  [logistic regression] calculation failed: %s" % e)

# ---------------------------------------------------------------
# [Method 1-B scaffold] bulk myeloid index partial correlation -- requires additional network fetch of myeloid markers
#   separately fetch CD68/CD14/LYZ/CSF1R/ITGAM into expr in fetch_myeloid_markers.py
#    Define functions here for the main process to call  <-- TRANSLATE
# ---------------------------------------------------------------
def main_B(myeloid_idx_df):
    """myeloid_idx_df: index=patient, columns=['myeloid_idx']; must include ZP3/TREM2 with same index"""
    log("=" * 70)
    log("[Check B] whether ZP3↔TREM2 at bulk level is independent of total myeloid burden")
    df = myeloid_idx_df.dropna()
    if len(df) < 30:
        log("  too few samples (%d), skip" % len(df)); return
    # 1) ZP3 vs myeloid index
    r1, p1 = stats.pearsonr(df['ZP3'], df['myeloid_idx'])
    log("  1) ZP3 vs myeloid index: r=%.3f p=%.3g" % (r1, p1))
    # 2) TREM2 vs myeloid index
    r2, p2 = stats.pearsonr(df['TREM2'], df['myeloid_idx'])
    log("  2) TREM2 vs myeloid index: r=%.3f p=%.3g" % (r2, p2))
    # 3) ZP3 vs TREM2 (uncorrected)
    r0, p0 = stats.pearsonr(df['ZP3'], df['TREM2'])
    log("  3) ZP3 vs TREM2 (crude correlation): r=%.3f p=%.3g" % (r0, p0))
    # 4) Partial correlation: controlling for myeloid index
    def partial_corr(x, y, cov):
        def reg(a, b):
            m = np.column_stack([np.ones_like(b), b])
            coef, *_ = np.linalg.lstsq(m, a, rcond=None)
            return coef[0] + m[:,1]*coef[1]
        rx = x - reg(x, cov)
        ry = y - reg(y, cov)
        r, p = stats.pearsonr(rx, ry)
        return r, p
    rp, pp = partial_corr(df['ZP3'].values, df['TREM2'].values, df['myeloid_idx'].values)
    log("  4) ZP3 vs TREM2 (partial correlation, controlling for myeloid index): r_partial=%.3f p=%.3g" % (rp, pp))
    if pp < 0.05:
        log("      Interpretation: after controlling for myeloid burden, still significant -> ZP3↔TREM2 has a specific association independent of total myeloid abundance")
    else:
        log("      Interpretation: after controlling for myeloid burden, not significant -> the apparent bulk ZP3↔TREM2 correlation is mainly driven by total myeloid burden (confounding)")
    r_pct = 100*(1 - pp/ (p0+1e-300))
    log("     (crude p=%.3g -> partial p=%.3g, if not significant it indicates that part is explained by total myeloid abundance)" % (p0, pp))

if __name__ == "__main__":
    main_A()
    print()
    main_B(None)  # silent when None (real data is passed in by the main script)
