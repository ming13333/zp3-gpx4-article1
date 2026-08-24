# -*- coding: utf-8 -*-
"""
GSE84465 official annotation analysis (2026-08-10 v3 correction)

Background: the v1/v2 GSE84465 analysis used marker scoring gating (pan-myeloid mean>0), classifying 62.4%
  cells as “myeloid”. After cross-checking with the official GSE series matrix annotation, this gating misclassified 302 official
  Neoplastic (tumor) cells + OPC/Vascular etc. into the myeloid pool, causing the intra-myeloid ZP3↔TREM2
  co-enrichment to appear as a false “reverse” (OR=0.57). This script recalculates with the official annotation and reaches a corrected conclusion.

Outputs:
  - gse84465_official_annotation.csv   （per-cell official type）
  - gse84465_official_coenrichment.csv （co-enrichment within the official pool）
  - gse84465_official_source.csv       （ZP3 cell-type enrichment）
  - Updated cross-cohort comparison with GSE182109
"""
import os, gzip, re
import numpy as np
import pandas as pd
import scipy.stats as stats

OUT = os.path.dirname(os.path.abspath(__file__))
SERIES_MATRIX = os.path.join(OUT, "gse84465_series_matrix.txt.gz")
EXPR_PATH = os.path.join(OUT, "..", "h1_pilot", "GSE84465_GBM_All_data.csv.gz")

def load_official_annotation():
    lines = []
    with gzip.open(SERIES_MATRIX, "rt") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!"):
                lines.append(line.rstrip())
    descs, ct_col = None, None
    for line in lines:
        m = re.match(r"!Sample_description\s+(.*)", line, re.S)
        if m:
            descs = [c.strip('"') for c in m.group(1).split("\t")]
        m = re.match(r"!Sample_characteristics_ch1\s+(.*)", line, re.S)
        if m:
            cols = [c.strip('"') for c in m.group(1).split("\t")]
            if any("cell type" in c for c in cols):
                ct_col = [c.replace("cell type: ", "") for c in cols]
    assert descs and ct_col and len(descs) == len(ct_col)
    return pd.DataFrame({"cell_id": descs, "official_type": ct_col})

def load_expr():
    with gzip.open(EXPR_PATH, "rt") as f:
        first = f.readline()
    sep = r"\s+" if first.count(" ") > first.count(",") else ","
    df = pd.read_csv(EXPR_PATH, sep=sep, index_col=0, compression="gzip")
    expr = df.T.apply(pd.to_numeric, errors="coerce").fillna(0)
    expr.index = [str(i) for i in expr.index]
    return expr

def fisher(a, b):
    a = np.asarray(a, bool); b = np.asarray(b, bool)
    aa = int((a & b).sum()); ab = int((a & ~b).sum())
    ba = int((~a & b).sum()); bb = int((~a & ~b).sum())
    OR, p = stats.fisher_exact([[aa, ab], [ba, bb]])
    return aa, ab, ba, bb, OR, p

def main():
    annot = load_official_annotation().set_index("cell_id")
    expr = load_expr()
    annot = annot.loc[expr.index]
    assert len(annot) == len(expr) and annot["official_type"].notna().all()
    annot.to_csv(os.path.join(OUT, "gse84465_official_annotation.csv"))

    zp3 = expr["ZP3"].values.astype(float)
    trem2 = expr["TREM2"].values.astype(float)
    ot = annot["official_type"].values.astype(str)
    zpos = zp3 > 0; tpos = trem2 > 0
    n = len(zp3)

    print("=" * 74)
    print("GSE84465 official annotation version: ZP3 ↔ TREM2 co-enrichment + source analysis (v3 correction)")
    print("  Official type composition:")
    vc = pd.Series(ot).value_counts()
    for t, c in vc.items():
        print("    %-15s %5d (%4.1f%%)" % (t, c, 100 * c / n))

    # Official Immune cell pool (total immune cell pool in GEO annotation)
    rows = []
    print("\n[Co-enrichment] Whether ZP3+ is enriched for TREM2+ (official annotation pool)")
    pools = [("All cells", np.ones(n, bool)),
             ("Official Immune cell", ot == "Immune cell"),
             ("Official Neoplastic", ot == "Neoplastic")]
    for label, mask in pools:
        aa, ab, ba, bb, OR, p = fisher(zpos[mask], tpos[mask])
        frac = aa / (aa + ab) if aa + ab else float("nan")
        bg = tpos[mask].mean()
        print("  [%s] n=%d | ZP3+&TREM2+=%d/%d | TREM2+ in ZP3+=%.1f%% (background %.1f%%) | OR=%.2f p=%.3g"
              % (label, int(mask.sum()), aa, aa + ab, 100 * frac, 100 * bg, OR, p))
        rows.append({"pool": label, "n": int(mask.sum()), "n_zp3pos": aa + ab,
                     "n_zp3pos_trem2pos": aa,
                     "frac_zp3pos_trem2pos": (round(frac, 4) if frac == frac else None),
                     "bg_trem2pos": round(float(bg), 4), "OR": round(OR, 2), "p": p})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "gse84465_official_coenrichment.csv"), index=False)

    # ZP3 source: enrichment per official cell type
    print("\n[Source] Whether ZP3+ cells are biased toward a particular official cell type")
    src = []
    for t in vc.index:
        m = ot == t
        oo = int((zpos & m).sum()); oo2 = int((zpos & ~m).sum())
        nn = int((~zpos & m).sum()); nn2 = int((~zpos & ~m).sum())
        if oo + oo2 > 0 and nn + nn2 > 0:
            OR, p = stats.fisher_exact([[oo, oo2], [nn, nn2]])
            frac = 100 * oo / max(1, oo + oo2)
            bg = 100 * m.mean()
            print("  %-15s: proportion in ZP3+=%.1f%% (background %.1f%%) | OR=%.2f p=%.3g"
                  % (t, frac, bg, OR, p))
            src.append({"cell_type": t, "frac_in_zp3pos": round(frac, 2),
                        "bg_pct": round(bg, 2), "OR": round(OR, 2), "p": p})
    pd.DataFrame(src).to_csv(os.path.join(OUT, "gse84465_official_source.csv"), index=False)

    # ZP3+ official cell type distribution
    print("\n[Distribution] Official cell types of ZP3+ cells (n=%d):" % int(zpos.sum()))
    print(pd.Series(ot[zpos]).value_counts().to_string())
    print("=" * 74)

if __name__ == "__main__":
    main()
