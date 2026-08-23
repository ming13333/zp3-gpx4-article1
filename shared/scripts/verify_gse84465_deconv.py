# -*- coding: utf-8 -*-
"""
方法1·单细胞来源甄别 — GSE84465 独立队列细胞层面复现（跨队列稳健性检验）

背景：GSE182109 子集 h5ad（10x 3'，HVG 空间）在 verify_sources_deconv.py 检验A
  显示 ZP3+ 髓系中 89.3% 共表达 TREM2+（OR=20.5, p=8.5e-11）。
本脚本在第二个独立 GBM 单细胞队列 GSE84465（Darmanis 2017, SMART-seq2 全长，
3589 细胞, genes×cells, raw count 尺度）重做同款检验，评估跨队列稳健性。

产出：
  - 共富集表：全细胞 / 髓系内 / MG-TAM-DC 亚群（Fisher exact OR+p）
  - 逻辑回归：ZP3 表达量 -> TREM2+（per-unit OR + 似然比 p）
  - ZP3 来源谱系分析：ZP3+ 细胞是否富集髓系？ZP3 高表达细胞的谱系 marker 特征
  - 与 GSE182109 主队列逐项对比表

2026-08-10  Craft 模式执行（补跨队列稳健性）
"""
import os, gzip
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.optimize import minimize

OUT = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(os.path.dirname(OUT), "h1_pilot", "GSE84465_GBM_All_data.csv.gz")
LOG = []

def log(msg=""):
    print(msg)
    LOG.append(str(msg))

# ---- marker 门控（与 h1_replicate_gse84465.py 完全一致，保证可比）----
pan_myeloid = ["CD68", "LYZ", "C1QA", "C1QB", "ITGAM", "CSF1R", "CD14"]
MG  = ["CX3CR1", "P2RY12", "TMEM119", "SALL1", "SIGLEC11"]
TAM = ["CD163", "VSIG4", "MRC1", "MSR1", "FOLR2"]
DC  = ["CLEC9A", "FCER1A", "CD1C", "LAMP3", "BATF3", "ITGAX"]

def fisher_table(a, b):
    a = np.asarray(a, bool); b = np.asarray(b, bool)
    aa = int(((a) & (b)).sum()); ab = int(((a) & (~b)).sum())
    ba = int(((~a) & (b)).sum()); bb = int(((~a) & (~b)).sum())
    OR, p = stats.fisher_exact([[aa, ab], [ba, bb]])
    return aa, ab, ba, bb, OR, p

def logistic_or(x, y):
    """logit(y)=b0+b1*x -> (OR=exp(b1), b0, b1)"""
    x = np.asarray(x, float); y = np.asarray(y, float)
    Xd = np.column_stack([np.ones_like(x), x])
    def negll(b):
        z = Xd @ b
        p = 1/(1+np.exp(-z))
        p = np.clip(p, 1e-12, 1-1e-12)
        return -(y*np.log(p) + (1-y)*np.log(1-p)).sum()
    res = minimize(negll, np.zeros(2), method='BFGS')
    return float(np.exp(res.x[1])), float(res.x[0]), float(res.x[1])

def lrt_p(x, y, b0, b1):
    """似然比检验 p（H0: b1=0）"""
    x = np.asarray(x, float); y = np.asarray(y, float)
    Xd = np.column_stack([np.ones_like(x), x])
    def ll(b):
        z = Xd @ b
        p = 1/(1+np.exp(-z))
        p = np.clip(p, 1e-12, 1-1e-12)
        return (y*np.log(p) + (1-y)*np.log(1-p)).sum()
    ll1 = ll([b0, b1])
    ll0 = ll([np.log(y.mean()/(1-y.mean())), 0.0])
    return float(stats.chi2.sf(2*(ll1 - ll0), 1))

def main():
    with gzip.open(PATH, "rt") as f:
        first = f.readline()
    sep = r"\s+" if first.count(" ") > first.count(",") else ","
    df = pd.read_csv(PATH, sep=sep, index_col=0, compression="gzip")
    expr = df.T.apply(pd.to_numeric, errors="coerce").fillna(0)

    log("=" * 74)
    log("GSE84465 细胞层面 ZP3↔TREM2 共富集 + ZP3 来源谱系（独立队列复现）")
    log("  cells=%d genes=%d | 平台: SMART-seq2 全长 | index=plate 位置（无官方注释）"
        % expr.shape)
    zp3 = expr["ZP3"].values.astype(float)
    trem2 = expr["TREM2"].values.astype(float)
    zp3_pos = zp3 > 0; trem2_pos = trem2 > 0
    n = len(zp3)
    log("  背景: ZP3+ =%d (%.2f%%) | TREM2+ =%d (%.2f%%)"
        % (zp3_pos.sum(), 100*zp3_pos.mean(), trem2_pos.sum(), 100*trem2_pos.mean()))
    log("  数据尺度检查: GAPDH p50=%.0f (raw count 量级) | ZP3 max=%.0f | TREM2 max=%.0f"
        % (np.median(expr["GAPDH"].values.astype(float)), zp3.max(), trem2.max()))

    # ---- marker 门控 ----
    def cmean(gl):
        ix = [g for g in gl if g in expr.columns]
        return expr[ix].mean(axis=1) if ix else pd.Series(0.0, index=expr.index)
    pm = cmean(pan_myeloid).values.astype(float)
    mg = cmean(MG); tam = cmean(TAM); dc = cmean(DC)
    myeloid = pm > 0
    sub = pd.DataFrame({"MG": mg, "TAM": tam, "DC": dc}).loc[myeloid]
    best = sub.idxmax(axis=1); val = sub.max(axis=1)
    subclass_pool = best.where(val > 0, "Unassigned")
    subclass_full = pd.Series("Not_myeloid", index=expr.index)
    subclass_full.loc[myeloid] = subclass_pool
    sc = subclass_full.values.astype(str)
    log("\n门控（与 GSE84465 H1 复现一致）: 泛髓系 n=%d (%.1f%%) | TAM=%d MG=%d DC=%d Un=%d"
        % (myeloid.sum(), 100*myeloid.mean(),
           int((sc == "TAM").sum()), int((sc == "MG").sum()),
           int((sc == "DC").sum()), int((sc == "Unassigned").sum())))

    # ---- 表1: 共富集 ----
    log("\n[表1] ZP3+ 是否富集 TREM2+")
    rows = []
    for label, mask in [("全细胞", np.ones(n, bool)), ("髓系内", myeloid)]:
        aa, ab, ba, bb, OR, p = fisher_table(zp3_pos[mask], trem2_pos[mask])
        frac = aa/(aa+ab) if aa+ab else float('nan')
        bg = (trem2_pos[mask]).mean()
        log("  [%s] n=%d | ZP3+&TREM2+=%d  ZP3+&TREM2-=%d | ZP3+ 中 TREM2+ 比例=%.1f%% (背景 %.1f%%) | OR=%.2f p=%.3g"
            % (label, int(mask.sum()), aa, ab, 100*frac, 100*bg, OR, p))
        rows.append({"level": label, "n": int(mask.sum()), "n_zp3pos": aa+ab,
                     "n_zp3pos_trem2pos": aa, "frac_zp3pos_trem2pos": round(frac, 4),
                     "bg_trem2pos": round(float(bg), 4), "OR": round(OR, 2), "p": p})
    for s in ["MG", "TAM", "DC"]:
        mask = sc == s
        if int(mask.sum()) < 10:
            continue
        aa, ab, ba, bb, OR, p = fisher_table(zp3_pos[mask], trem2_pos[mask])
        frac = aa/(aa+ab) if aa+ab else float('nan')
        bg = (trem2_pos[mask]).mean()
        log("  [%s] n=%d | ZP3+&TREM2+=%d/%d | 亚群内 TREM2+ 背景=%.1f%% | OR=%.2f p=%.3g"
            % (s, int(mask.sum()), aa, aa+ab, 100*bg, OR, p))
        rows.append({"level": s, "n": int(mask.sum()), "n_zp3pos": aa+ab,
                     "n_zp3pos_trem2pos": aa,
                     "frac_zp3pos_trem2pos": (round(frac, 4) if frac == frac else None),
                     "bg_trem2pos": round(float(bg), 4), "OR": round(OR, 2), "p": p})

    # ---- 表2: 逻辑回归 ----
    log("\n[表2] ZP3 表达量 -> TREM2+（log1p, 似然比检验）")
    lr_rows = []
    x_all = np.log1p(zp3); y_all = trem2_pos.astype(float)
    or_a, b0, b1 = logistic_or(x_all, y_all)
    p_a = lrt_p(x_all, y_all, b0, b1)
    log("  全细胞: OR=%.2f / 单位 log1p ZP3, p=%.3g" % (or_a, p_a))
    lr_rows.append({"level": "全细胞", "per_unit_OR": round(or_a, 2), "p": p_a})
    if myeloid.sum() > 0:
        x_m = np.log1p(zp3[myeloid]); y_m = trem2_pos[myeloid].astype(float)
        or_m, b0m, b1m = logistic_or(x_m, y_m)
        p_m = lrt_p(x_m, y_m, b0m, b1m)
        log("  髓系内: OR=%.2f / 单位 log1p ZP3, p=%.3g" % (or_m, p_m))
        lr_rows.append({"level": "髓系内", "per_unit_OR": round(or_m, 2), "p": p_m})

    # ---- 表3: ZP3 是否髓系来源（独立于 TREM2 的直接检验）----
    log("\n[表3] ZP3+ 细胞是否富集髓系（'bulk ZP3 来自髓系?'的直接检验）")
    src_rows = []
    for pth in [0, 5]:
        a = zp3 > pth
        if int(a.sum()) < 10:
            continue
        aa, ab, ba, bb, OR, p = fisher_table(a, myeloid)
        frac = aa/(aa+ab) if aa+ab else float('nan')
        bg = myeloid.mean()
        log("  ZP3>%d: ZP3+ 中髓系=%.1f%% (背景 %.1f%%) | OR=%.2f p=%.3g (n=%d)"
            % (pth, 100*frac, 100*bg, OR, p, int(a.sum())))
        src_rows.append({"zp3_threshold": pth, "n_zp3pos": int(a.sum()),
                         "frac_zp3pos_myeloid": round(frac, 4),
                         "bg_myeloid": round(float(bg), 4),
                         "OR": round(OR, 2), "p": p})

    # ---- 表4: ZP3 高表达细胞谱系 marker 特征 ----
    log("\n[表4] ZP3 高表达细胞 (ZP3>10, n=%d) vs ZP3=0 (n=%d) 谱系 marker log2FC"
        % (int((zp3 > 10).sum()), int((zp3 == 0).sum())))
    markers = {
        "肿瘤/神经": ["EGFR", "OLIG2", "SOX2", "GFAP", "NES", "TOP2A"],
        "髓系": ["CD68", "LYZ", "CD163", "CX3CR1", "P2RY12", "TMEM119"],
        "T/NK": ["NKG7", "CD3D"],
        "内皮/周皮": ["VWF", "PDGFRB", "COL1A1"],
    }
    hi_mask = (zp3 > 10); lo_mask = (zp3 == 0)
    marker_rows = []
    for grp, gl in markers.items():
        avail = [g for g in gl if g in expr.columns]
        if not avail:
            continue
        hi = expr.loc[hi_mask, avail].mean(axis=0)
        lo = expr.loc[lo_mask, avail].mean(axis=0)
        lfc = np.log2((hi.values + 1) / (lo.values + 1))
        for i, g in enumerate(avail):
            log("  %-8s (%-6s): hi=%.1f lo=%.1f log2FC=%+.2f"
                % (g, grp, hi.iloc[i], lo.iloc[i], lfc[i]))
            marker_rows.append({"gene": g, "group": grp,
                                "zp3hi_mean": round(float(hi.iloc[i]), 2),
                                "zp3lo_mean": round(float(lo.iloc[i]), 2),
                                "log2FC": round(float(lfc[i]), 2)})

    # ---- 保存 ----
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "gse84465_coenrichment.csv"), index=False)
    pd.DataFrame(lr_rows).to_csv(os.path.join(OUT, "gse84465_logistic.csv"), index=False)
    pd.DataFrame(src_rows).to_csv(os.path.join(OUT, "gse84465_source_myeloid.csv"), index=False)
    pd.DataFrame(marker_rows).to_csv(os.path.join(OUT, "gse84465_zp3hi_markers.csv"), index=False)
    log("\n已保存: gse84465_coenrichment.csv / gse84465_logistic.csv / gse84465_source_myeloid.csv / gse84465_zp3hi_markers.csv")
    log("=" * 74)

if __name__ == "__main__":
    main()
