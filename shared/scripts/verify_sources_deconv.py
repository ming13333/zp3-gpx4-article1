# -*- coding: utf-8 -*-
"""
方法1（单细胞反卷积/来源甄别）+ 响应"bulk ZP3 来源 = 肿瘤 vs 髓系混杂"的审稿质疑。

核心思路（不依赖脆弱的细胞类型参考谱反卷积，改用两层针对性检验，直击混杂质疑）：

  [检验 A] 单细胞层面 ZP3↔TREM2 共富集（本地 h5ad，零网络）
     - 直接检验：ZP3+ 细胞是否富集 TREM2+（全细胞 / 髓系内 / MG-TAM-DC 亚群内）
     - 若 ZP3+ 在髓系内几乎全部共表达 TREM2，则"bulk ZP3↔TREM2 相关"不是
       肿瘤细胞混杂所致，而是同一群 TREM2+ 髓系细胞同时贡献 ZP3 与 TREM2。
     - 用 Fisher 精确检验 + 富集强度（OR）+ 单变量逻辑回归（logistic OR per unit ZP3）

  [检验 B] bulk 层面：ZP3 与 TREM2 的关联是否独立于"总髓系负荷"
     - 用泛髓系 marker（CD68/CD14/LYZ/CSF1R/ITGAM）构建髓系指数
     - Pearson: ZP3 vs 髓系指数
     - 偏相关: 控制髓系指数后 ZP3 vs TREM2（用残差相关实现）
     - 若控制髓系负荷后关联大幅衰减 -> 混杂（bulk 层面 ZP3↔TREM2 表观相关来自髓系总量）
     - 若衰减但仍显著 -> 存在独立于总髓系量的特异关联
     - 两种结果都是"可发表"的诚实结论，但解读不同。

  2026-08-10  Craft 模式执行
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
# [检验 A] 单细胞层面 ZP3 ↔ TREM2 共富集
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
    log("[检验 A] 单细胞层面 ZP3+ 细胞是否富集 TREM2+  (n_cells=%d)" % len(zp3))
    log("  var 空间 = %d HVG（非全基因组）；ZP3/TREM2 均因高变异被保留" % adata.n_vars)
    log("  背景: ZP3+ %d 细胞 (%.2f%%), TREM2+ %d 细胞 (%.2f%%)"
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
        log("      ZP3+ 中 TREM2+ 比例=%.1f%% | 整体背景 TREM2+ 比例=%.1f%%" % (100*frac_zp3_trem2, 100*frac_bg_trem2))
        log("      Fisher OR=%.2f, p=%.3g" % (OR, p))
        return OR, p

    # A0 全细胞
    fisher_or_table(zp3_pos, trem2_pos, "全细胞")
    # A1 髓系内
    myeloid = obs['myeloid'].values.astype(bool)
    if myeloid.sum() > 0:
        m = fisher_or_table(zp3_pos[myeloid], trem2_pos[myeloid], "髓系内 myel=True")
    # A2 按髓系亚群分层
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
                log("  [亚群:%s n=%d]  ZP3+&TREM2+=%d/%d  | 亚群内 TREM2+ 背景=%.1f%% | Fisher OR=%.2f p=%.3g"
                    % (s, sel.sum(), aa, aa+ab, 100*ftrem.mean(), OR, p))
    # A3 单变量逻辑回归（ZP3 表达量 -> TREM2+，全细胞）给出每单位 ZP3 的 OR
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
        log("  [逻辑回归·未校正] logit(TREM2+)=%.3f+%.3f*ZP3 -> ZP3 每 +1(log1p) 的 OR=%.2f"
            % (b0, b1, np.exp(b1)))
    except Exception as e:
        log("  [逻辑回归] 计算失败: %s" % e)

# ---------------------------------------------------------------
# [方法1-B 支架] bulk 髓系指数偏相关 —— 需要额外网络拉取髓系 marker
#   单独在 fetch_myeloid_markers.py 中拉取 CD68/CD14/LYZ/CSF1R/ITGAM 到 expr
#   此处定义函数供主流程调用
# ---------------------------------------------------------------
def main_B(myeloid_idx_df):
    """myeloid_idx_df: index=patient, 列=['myeloid_idx']；需含 ZP3/TREM2 同 index"""
    log("=" * 70)
    log("[检验 B] bulk 层面 ZP3↔TREM2 是否独立于总髓系负荷")
    df = myeloid_idx_df.dropna()
    if len(df) < 30:
        log("  样本太少(%d)，跳过" % len(df)); return
    # 1) ZP3 vs 髓系指数
    r1, p1 = stats.pearsonr(df['ZP3'], df['myeloid_idx'])
    log("  1) ZP3 vs 髓系指数: r=%.3f p=%.3g" % (r1, p1))
    # 2) TREM2 vs 髓系指数
    r2, p2 = stats.pearsonr(df['TREM2'], df['myeloid_idx'])
    log("  2) TREM2 vs 髓系指数: r=%.3f p=%.3g" % (r2, p2))
    # 3) ZP3 vs TREM2 (未校正)
    r0, p0 = stats.pearsonr(df['ZP3'], df['TREM2'])
    log("  3) ZP3 vs TREM2 (粗相关): r=%.3f p=%.3g" % (r0, p0))
    # 4) 偏相关：控制髓系指数
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
    log("  4) ZP3 vs TREM2 (偏相关，控制髓系指数): r_partial=%.3f p=%.3g" % (rp, pp))
    if pp < 0.05:
        log("     解读:控制髓系负荷后仍显著 -> ZP3↔TREM2 存在独立于总髓系量的特异关联")
    else:
        log("     解读:控制髓系负荷后不显著 -> bulk ZP3↔TREM2 表观相关主要由总髓系负荷驱动(混杂)")
    r_pct = 100*(1 - pp/ (p0+1e-300))
    log("     (粗p=%.3g -> 偏p=%.3g, 若不显著则说明部分由髓系总量解释)" % (p0, pp))

if __name__ == "__main__":
    main_A()
    print()
    main_B(None)  # None 时静默（真实数据由主脚本传入）
