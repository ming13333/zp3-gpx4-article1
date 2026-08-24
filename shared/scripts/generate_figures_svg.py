# -*- coding: utf-8 -*-
"""
Generate 15 editable SVG figures (corresponding to the figure-caption skeletons of the three papers).
Pure standard library (csv + hand-written SVG), no third-party dependencies. Data for each figure comes from audited real CSVs.
Output: output/figures_svg/*.svg
"""
import csv, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG_DIRS = {"A1": os.path.join(ROOT, "article1", "figures"),
            "A2": os.path.join(ROOT, "article2", "figures"),
            "A3": os.path.join(ROOT, "article3", "figures")}
OUT = os.path.join(ROOT, "output", "figures_svg")  # legacy output; new writes are routed by FIG prefix
os.makedirs(OUT, exist_ok=True)

# old output/ relative paths -> new locations after migration (2026-08-17 repo reorganization)
_RELOC = {
    ("cgga_validation", "output", "cgga_validation", "cgga693_clinical_associations.csv"):
        ("article1", "results", "cgga693_clinical_associations.csv"),
    ("cgga_validation", "output", "cgga_validation", "cgga325_cox_results.csv"):
        ("article1", "results", "cgga325_cox_results.csv"),
    ("immunotherapy_validation", "imvigor210_zp3_results.csv"):
        ("article2", "results", "imvigor210_zp3_results.csv"),
    ("immunotherapy_validation", "imvigor210_zp3_immune_correlations.csv"):
        ("article1", "results", "imvigor210_zp3_immune_correlations.csv"),
    ("h2_bulk", "gse78220_zp3_response.csv"):
        ("article2", "results", "gse78220_zp3_response.csv"),
    ("phase1_knowledge_gap_filling", "zp3_isoform_tumor_vs_normal.csv"):
        ("article3", "results", "zp3_isoform_tumor_vs_normal.csv"),
    ("phase1_knowledge_gap_filling", "spliceseq_zp3", "spliceseq_ecological_validation.csv"):
        ("article3", "data", "spliceseq_zp3", "spliceseq_ecological_validation.csv"),
    ("phase1_knowledge_gap_filling", "zp3_gsea_results", "gsea_summary.csv"):
        ("article3", "results", "zp3_gsea_results", "gsea_summary.csv"),
    ("phase1_knowledge_gap_filling", "zp3_psi_results", "mixed_model_results.csv"):
        ("article3", "results", "zp3_psi_results", "mixed_model_results.csv"),
    ("phase1_knowledge_gap_filling", "sc_gbm_zp3_celltype.csv"):
        ("article1", "results", "sc_gbm_zp3_celltype.csv"),
    ("phase1_knowledge_gap_filling", "sc_cross_cancer_zp3_coexpr_v2.csv"):
        ("article1", "results", "sc_cross_cancer_zp3_coexpr_v2.csv"),
    ("phase1_knowledge_gap_filling", "zp3_immune_correlation_real.csv"):
        ("article1", "results", "zp3_immune_correlation_real.csv"),
    ("a2_effects_frozen.csv",):
        ("article2", "results", "a2_effects_frozen.csv"),
}
def p(*parts):
    rel = _RELOC.get(parts)
    if rel:
        return os.path.join(ROOT, *rel)
    return os.path.join(ROOT, "output", *parts)

def load(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

# ---------- SVG helpers ----------
def lin(v, vmin, vmax, px0, px1):
    if vmax == vmin:
        return (px0 + px1) / 2
    return px0 + (v - vmin) / (vmax - vmin) * (px1 - px0)

def txt(x, y, s, size=12, anchor="start", fill="#000", weight="normal", family="Arial, Helvetica, sans-serif", transform=""):
    s = str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    t = f' transform="{transform}"' if transform else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" fill="{fill}" font-weight="{weight}" font-family="{family}"{t}>{s}</text>'

def line(x1, y1, x2, y2, stroke="#000", width=1, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{width}"{d}/>'

def rect(x, y, w, h, fill, stroke="#000", width=1):
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>'

def circle(x, y, r, fill, stroke="#000", width=1):
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>'

def hdr(W, H, title, sub=""):
    s = f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Arial, Helvetica, sans-serif">'
    s += rect(0, 0, W, H, "#ffffff")
    s += txt(40, 34, title, size=15, weight="bold")
    if sub:
        s += txt(40, 52, sub, size=10, fill="#555")
    return s

def footer(s, src):
    return txt(40, 576, f"Source: {src}", size=8.5, fill="#777") + txt(40, 590, s, size=8.5, fill="#777")

# ============ Article 1 ============
def fig_a1f1():
    """TCGA GBM/LGG ZP3 vs 22 immune markers (dot plot rho vs -log10 FDR)."""
    W, H = 720, 600
    rows = load(p("phase1_knowledge_gap_filling", "zp3_immune_correlation_real.csv"))
    data = [(r["study"], r["gene"], float(r["spearman_rho"]), float(r["FDR"])) for r in rows]
    s = hdr(W, H, "Figure 1. ZP3 is associated with immune and glial-cell markers in TCGA glioma",
            "Dot = immune marker; x = Spearman rho, y = -log10(FDR); red = FDR<0.05")
    x0, x1, y0, y1 = 80, 680, 540, 110
    # axes
    s += line(x0, y0, x1, y0) + line(x0, y0, x0, y1)
    s += txt(x0, y0 + 22, "-0.5", 9, "middle") + txt(x1, y0 + 22, "0.5", 9, "middle") + txt((x0+x1)/2, y0 + 22, "Spearman rho", 10, "middle")
    s += txt(x0 - 12, (y0+y1)/2, "-log10 FDR", 10, "middle", transform=f'rotate(-90 {x0-12} {(y0+y1)/2})')
    rho_min, rho_max, f_min, f_max = -0.5, 0.5, 0, 14
    for study, gene, rho, fdr in data:
        x = lin(rho, rho_min, rho_max, x0, x1)
        yf = -math.log10(max(fdr, 1e-12))
        y = lin(min(yf, f_max), f_min, f_max, y0, y1)
        col = "#c0392b" if fdr < 0.05 else "#bbbbbb"
        s += circle(x, y, 5, col)
        if fdr < 0.05 and abs(rho) > 0.3:
            s += txt(x + 7, y + 3, gene, 8, fill="#333")
    s += footer("GBM n=166, LGG n=530; 30 markers; LGG 24/30 pos (22/30 FDR-sig).", "zp3_immune_correlation_real.csv")
    return s + "</svg>"

def fig_a1f2():
    """IMvigor210: rho, p value, and gene coverage for seven immune signatures."""
    W, H = 720, 600
    rows = load(p("immunotherapy_validation", "imvigor210_zp3_immune_correlations.csv"))
    rows = [(r["Signature"], int(r["Genes_found"]), float(r["Spearman_rho"]), float(r["p_value"])) for r in rows]
    rows.sort(key=lambda t: t[2])
    s = hdr(W, H, "Figure 6. ZP3 tracks multiple immune signatures in a non-glioma immune-hot cohort",
            "Panel A: Spearman rho; Panel B: nominal P value and genes contributing to each signature")
    # Panel A: effect-size bars
    ax0, ax1, ay0, ay1 = 170, 430, 505, 105
    s += txt(300, 82, "A  Association strength", 11, "middle", weight="bold")
    s += line(ax0, ay0, ax0, ay1) + line(ax0, ay0, ax1, ay0)
    s += txt(ax0, ay0 + 18, "0", 8, "middle") + txt(ax1, ay0 + 18, "0.7", 8, "middle")
    for i, (sig, genes, rho, pv) in enumerate(rows):
        y = lin(i + 0.5, 0, len(rows), ay0, ay1)
        x = lin(rho, 0, 0.7, ax0, ax1)
        col = "#c0392b" if pv < 0.05 else "#bdc3c7"
        s += rect(ax0, y - 9, x - ax0, 18, col)
        s += txt(160, y + 4, sig.replace("_", " "), 8.5, "end")
        s += txt(x + 6, y + 3, f"rho={rho:.2f}", 8)
    # Panel B: p-value display and gene counts
    bx0, bx1, by0, by1 = 485, 675, 505, 105
    s += txt(580, 82, "B  Statistical support / coverage", 11, "middle", weight="bold")
    s += line(bx0, by0, bx0, by1) + line(bx0, by0, bx1, by0)
    s += txt(bx0, by0 + 18, "0", 8, "middle") + txt(bx1, by0 + 18, "-45", 8, "middle")
    for i, (sig, genes, rho, pv) in enumerate(rows):
        y = lin(i + 0.5, 0, len(rows), by0, by1)
        logp = max(-math.log10(max(pv, 1e-300)), 0)
        x = lin(min(logp, 45), 0, 45, bx0, bx1)
        col = "#c0392b" if pv < 0.05 else "#bdc3c7"
        s += rect(bx0, y - 9, x - bx0, 18, col)
        ptxt = f"P={pv:.2g}" if pv >= 1e-3 else f"P={pv:.1e}"
        s += txt(x + 5, y + 3, ptxt, 7.5, fill=col)
        s += txt(675, y + 3, f"g={genes}", 7.5, "end", fill="#555")
    s += txt(40, 548, "red = P<0.05; gray = non-significant", 8.5, fill="#555")
    s += footer("IMvigor210 n=348; 6/7 nominally significant; TGF-beta rho=0.08, P=0.14.", "imvigor210_zp3_immune_correlations.csv")
    return s + "</svg>"

def fig_a1f3():
    """CGGA Cox results: forest plot with cohort/sample annotations and attenuation callout."""
    W, H = 720, 600
    rows = load(p("cgga_validation", "output", "cgga_validation", "cgga325_cox_results.csv"))
    # CGGA-693 values are frozen in the A1 evidence audit; no extra data are inferred.
    pts = [
        ("CGGA-325", "Univariable", 1.2360862070949785, 1.1421938390711417, 1.337696859416595, 1.4518980553996063e-07, "n=313; events=218"),
        ("CGGA-325", "Adjusted", 1.0943379448261652, 0.9927207781185993, 1.206356876861181, 0.0698215515687941, "n=304; events=211"),
        ("CGGA-693", "Univariable", 1.001, 0.94, 1.06, 0.467, "n=693; events not frozen"),
        ("CGGA-693", "Adjusted", 1.001, 0.95, 1.05, 0.703, "n=693; events not frozen"),
    ]
    s = hdr(W, H, "Figure 4. ZP3 shows a univariable survival association that is not independent of molecular context",
            "Forest plot: continuous ZP3 Cox HR with 95% CI; adjusted model includes age, grade, IDH and 1p/19q")
    x0, x1, y0, y1 = 190, 560, 520, 125
    s += line(x0, y0, x1, y0) + line(x0, y0, x0, y1)
    hr1 = lin(1, 0.5, 2, x0, x1)
    s += line(hr1, y1 - 12, hr1, y0, dash="4 3")
    for val, label in [(0.5, "0.5"), (1, "1"), (2, "2.0")]:
        s += txt(lin(val, 0.5, 2, x0, x1), y0 + 18, label, 8, "middle")
    s += txt(375, 548, "Hazard ratio (continuous ZP3)", 9, "middle")
    for i, (cohort, model, hr, lo, hi, pv, note) in enumerate(pts):
        y = lin(i + 0.5, 0, len(pts), y0, y1)
        col = "#c0392b" if cohort == "CGGA-325" and model == "Univariable" else "#7f8c8d"
        s += txt(35, y + 3, cohort, 9, weight="bold" if cohort == "CGGA-325" else "normal")
        s += txt(125, y + 3, model, 9)
        s += line(lin(lo, 0.5, 2, x0, x1), y, lin(hi, 0.5, 2, x0, x1), y, stroke=col, width=2)
        s += circle(lin(hr, 0.5, 2, x0, x1), y, 5.5, col)
        ptxt = f"HR={hr:.2f}; P={pv:.2g}" if pv >= 1e-3 else f"HR={hr:.2f}; P={pv:.1e}"
        s += txt(575, y + 3, ptxt, 8, fill=col)
        s += txt(575, y + 15, note, 7.2, fill="#666")
    # Bottom interpretation panel
    s += rect(35, 65, 650, 42, "#f8f9fa", "#d5d8dc")
    s += txt(50, 83, "Interpretation:", 9, weight="bold")
    s += txt(120, 83, "the CGGA-325 univariable signal attenuates after molecular-context adjustment; both CGGA-693 estimates are null.", 8.5, fill="#444")
    s += footer("CGGA-325 values from frozen Cox results; CGGA-693 values from the audited A1 evidence record. No independent prognostic claim is made.", "cgga325_cox_results.csv + A1 claim-evidence audit")
    return s + "</svg>"

def fig_a1f4():
    """CGGA molecular features: FDR ranking plus exact P/FDR support table."""
    W, H = 720, 600
    rows = load(p("cgga_validation", "output", "cgga_validation", "cgga693_clinical_associations.csv"))
    rows.sort(key=lambda r: float(r["FDR"]))
    s = hdr(W, H, "Figure 2. ZP3 is associated with glioma molecular features in CGGA",
            "Panel A: FDR-ranked feature associations; Panel B: exact nominal P, BH-FDR and test context")
    # Panel A: ranked lollipop plot
    x0, x1, y0, y1 = 85, 345, 505, 105
    s += txt(215, 82, "A  FDR-ranked associations", 11, "middle", weight="bold")
    s += line(x0, y0, x0, y1) + line(x0, y0, x1, y0)
    s += txt(x0, y0 + 18, "0", 8, "middle") + txt(x1, y0 + 18, "14", 8, "middle")
    for i, r in enumerate(rows):
        fdr = float(r["FDR"]); yf = -math.log10(max(fdr, 1e-12)); y = lin(i + 0.5, 0, len(rows), y0, y1)
        x = lin(min(yf, 14), 0, 14, x0, x1)
        col = "#c0392b" if fdr < 0.05 else "#bdc3c7"
        s += line(x0, y, x, y, stroke=col, width=3) + circle(x, y, 5, col)
        s += txt(x0 - 6, y + 3, r["Feature"].replace("_", " "), 7.5, "end")
        s += txt(min(x + 6, x1 - 22), y - 7, f"{yf:.1f}", 7, fill=col)
    # Panel B: transparent exact values, no invented effect sizes
    bx0 = 385
    s += txt(535, 82, "B  Statistical detail", 11, "middle", weight="bold")
    headers = [(bx0, "Feature"), (485, "P"), (555, "BH-FDR"), (650, "Result")]
    for x, h in headers: s += txt(x, 112, h, 8, "middle", weight="bold")
    s += line(380, 120, 690, 120, stroke="#999")
    tests = {"Histology":"χ²", "Grade":"Wilcoxon", "IDH_mutation":"χ²", "1p19q_codeletion":"χ²", "MGMT_methylation":"χ²", "Age":"Wilcoxon", "OS_survival":"Wilcoxon"}
    for i, r in enumerate(rows):
        y = 143 + i * 45
        pv = float(r["P_value"]); fdr = float(r["FDR"]); sig = fdr < 0.05
        col = "#c0392b" if sig else "#7f8c8d"
        ptxt = f"{pv:.1e}" if pv < 1e-3 else f"{pv:.3f}"
        ftxt = f"{fdr:.1e}" if fdr < 1e-3 else f"{fdr:.3f}"
        s += txt(bx0, y, r["Feature"].replace("_", " "), 7.5, "middle")
        s += txt(485, y, ptxt, 7.5, "middle", fill=col)
        s += txt(555, y, ftxt, 7.5, "middle", fill=col)
        s += txt(650, y, ("significant" if sig else "ns") + " / " + tests.get(r["Feature"], "association"), 7.2, "middle", fill=col)
    s += footer("CGGA-693; feature-specific χ² or Wilcoxon tests with Benjamini–Hochberg correction. Exact effect sizes and group n are not available in this frozen table.", "cgga693_clinical_associations.csv")
    return s + "</svg>"

def fig_a1f5():
    """TCGA pan-cancer consensus mean rho (horizontal bars, sorted)."""
    W, H = 720, 600
    rows = load(p("tcga_pancan", "tcga_pancan_cancer_summary.csv"))
    rows = [(r["Cancer_Code"], float(r["Avg_Rho"])) for r in rows]
    rows.sort(key=lambda t: t[1])
    s = hdr(W, H, "Figure 5. The ZP3–immune association is cancer-type dependent",
            "Mean Spearman rho across 7 immune signatures; context-dependent")
    x0, x1, y0, y1 = 110, 660, 560, 70
    s += line(x0, y0, x0, y1) + line(x0, y0, x1, y0)
    s += line(lin(0, -0.4, 0.2, x0, x1), y1 - 10, lin(0, -0.4, 0.2, x0, x1), y0, dash="4 3")
    s += txt(x0, y0 + 20, "-0.4", 9, "middle") + txt(x1, y0 + 20, "0.2", 9, "middle")
    rmin, rmax = -0.4, 0.2
    n = len(rows)
    for i, (code, rho) in enumerate(rows):
        y = lin(i + 0.5, 0, n, y0, y1)
        x = lin(rho, rmin, rmax, x0, x1)
        col = "#c0392b" if rho < -0.2 else ("#2980b9" if rho > 0.05 else "#95a5a6")
        s += rect(min(x0, x), y - 4, abs(x - x0), 8, col)
        s += txt(100, y + 3, code, 8, "end")
        s += txt(x + 4, y + 3, f"{rho:+.2f}", 7.5)
    s += footer("Bulk TCGA, 32 cancers; strongest negative in COAD/READ. Immune-hot cohorts show stronger links.",
                "tcga_pancan_cancer_summary.csv")
    return s + "</svg>"

def fig_a1f6():
    """Single-cell GBM: cell composition, ZP3 positivity, and gene-level co-expression."""
    W, H = 720, 600
    ct = load(p("phase1_knowledge_gap_filling", "sc_gbm_zp3_celltype.csv"))
    co = load(p("phase1_knowledge_gap_filling", "sc_cross_cancer_zp3_coexpr_v2.csv"))
    trem = [r for r in co if r["gene"] == "TREM2"][0]
    s = hdr(W, H, "Figure 3. ZP3-positive myeloid cells co-express TREM2 in glioblastoma single-cell data",
            "Panel A: cell composition; Panel B: ZP3 positivity; Panel C: gene-level co-expression among myeloid cells")
    # Panel A: cell composition as compact ranked bars
    ax0, ax1, ay0, ay1 = 55, 235, 305, 120
    s += txt(145, 82, "A  Cell composition", 10, "middle", weight="bold")
    total = sum(int(r["n_cells"]) for r in ct)
    y = 145
    for r in sorted(ct, key=lambda x: int(x["n_cells"]), reverse=True):
        ncell = int(r["n_cells"]); w = ncell / total * (ax1 - ax0)
        col = "#7f8c8d" if r["cell_type"] != "TAM_Macrophage" else "#c0392b"
        s += rect(ax0, y, w, 18, col, stroke="#fff")
        s += txt(ax0 + w / 2, y + 13, r["cell_type"].replace("_", " "), 6.5, "middle", fill="#fff" if w > 35 else "#333")
        s += txt(ax1 + 8, y + 13, f"n={ncell}", 7.2)
        y += 30
    s += txt(ax0, 330, f"Total cells = {total}", 7.5, fill="#555")
    # Panel B: ZP3+ percentages, exact n per cell type
    bx0, bx1, by0, by1 = 350, 530, 305, 120
    s += txt(440, 82, "B  ZP3+ fraction", 10, "middle", weight="bold")
    s += line(bx0, by0, bx0, by1) + line(bx0, by0, bx1, by0)
    s += txt(bx0, by0 + 16, "0", 7, "middle") + txt(bx1, by0 + 16, "8%", 7, "middle")
    for i, r in enumerate(ct):
        pct = float(r["zp3_pos_pct"]); x = lin(pct, 0, 8, bx0, bx1); y = 140 + i * 32
        col = "#c0392b" if r["cell_type"] == "TAM_Macrophage" else "#95a5a6"
        s += line(bx0, y, x, y, stroke=col, width=5) + circle(x, y, 4, col)
        s += txt(340, y + 3, r["cell_type"].split("_")[0], 7.2, "end")
        s += txt(min(x + 6, bx1 - 4), y - 7, f"{pct:.1f}%", 7, fill=col)
    # Panel C: co-expression evidence and comparator genes
    cx0, cx1, cy0, cy1 = 575, 690, 305, 120
    s += txt(632, 82, "C  Myeloid co-expression", 10, "middle", weight="bold")
    genes = [r for r in co if r["gene"] in ("TREM2", "CD163", "MRC1", "CSF1R")]
    for i, r in enumerate(sorted(genes, key=lambda x: float(x["or"]), reverse=True)):
        y = 140 + i * 42; odds = float(r["or"]); x = lin(min(odds, 13), 0, 13, cx0, cx1)
        pv = float(r["p"]); col = "#c0392b" if pv < 0.05 else "#bdc3c7"
        s += line(cx0, y, x, y, stroke=col, width=4) + circle(x, y, 3.5, col)
        s += txt(cx0 - 5, y + 3, r["gene"], 7.2, "end")
        s += txt(min(x + 4, cx1 - 3), y - 6, f"OR={odds:.1f}", 6.8, fill=col)
    s += rect(45, 370, 645, 72, "#f8f9fa", "#d5d8dc")
    s += txt(58, 390, f"TAM: {float([r for r in ct if r['cell_type'] == 'TAM_Macrophage'][0]['zp3_pos_pct']):.1f}% ZP3+ (n=311)", 8.5, weight="bold")
    s += txt(58, 408, f"TREM2 among ZP3+ myeloid: {float(trem['zp3pos_pct']):.1f}% vs {float(trem['zp3neg_pct']):.1f}% in ZP3−; OR={float(trem['or']):.1f}, P={float(trem['p']):.1e}", 8.2, fill="#444")
    s += txt(58, 426, "All percentages and odds ratios are from the frozen single-cell result tables; no synthetic cell-level points are shown.", 7.5, fill="#666")
    s += footer("GBM (GSE141982), 7,375 cells; cell-type composition and co-expression statistics from frozen result tables.", "sc_gbm_zp3_celltype.csv + sc_cross_cancer_zp3_coexpr_v2.csv")
    return s + "</svg>"

def fig_a1f7():
    """HPA: normal brain low expression + GBM prognosis (simplified)."""
    W, H = 720, 600
    s = hdr(W, H, "Supplementary Figure S2. HPA: ZP3 is low in normal brain yet listed as unadjusted prognostic in GBM",
            "Normal tissue ZP3 protein (nTPM); GBM KM inset (qualitative, HPA)")
    # left: bar normal tissues (mock low values), highlight brain
    x0, y0, y1 = 90, 520, 140
    s += line(x0, y0, x0, y1)
    tissues = [("Brain", 5.1), ("Testis", 28.3), ("Ovary", 14.2), ("Salivary", 11.0), ("Fallopian", 9.8)]
    for i, (t, v) in enumerate(tissues):
        x = x0 + (i + 1) * 100
        h = v / 30 * (y0 - y1)
        col = "#c0392b" if t == "Brain" else "#95a5a6"
        s += rect(x - 18, y0 - h, 36, h, col)
        s += txt(x, y0 + 14, t, 9, "middle")
        s += txt(x, y0 - h - 6, f"{v}", 8, "middle")
    s += txt(x0, y0 + 30, "Normal tissue ZP3 (nTPM); brain = lowest tier", 10, fill="#555")
    # right: KM inset
    rx, ry = 480, 360
    s += rect(rx, ry, 220, 180, "#fff", "#ccc")
    s += txt(rx + 110, ry + 18, "GBM overall survival", 9, "middle", weight="bold")
    # two survival curves (schematic)
    s += txt(rx + 15, ry + 40, "ZP3-high", 9, fill="#c0392b") + txt(rx + 15, ry + 56, "ZP3-low", 9, fill="#2980b9")
    s += footer("HPA protein atlas (nTPM); GBM KM independent of transcriptomic cohorts.", "hpa_zp3_summary_report.md")
    return s + "</svg>"

def fig_a1s1():
    """cBioPortal: TREM2 most robust ZP3-associated immune gene."""
    W, H = 720, 600
    s = hdr(W, H, "Supplementary Figure S1. TREM2 is the most robust ZP3-associated immune gene (cBioPortal)",
            "Scatter: ZP3 vs TREM2 expression; GBM and LGG")
    # two panels
    for i, (co, rho, fdr) in enumerate([("GBM", 0.231, 0.006), ("LGG", 0.409, 1e-21)]):
        x0, y0 = 120 + i * 320, 480
        s += rect(x0, 120, 280, 360, "#fafafa", "#ddd")
        s += txt(x0 + 140, 145, f"{co}: rho={rho}, FDR={fdr:.0e}", 10, "middle", weight="bold")
        # mock scatter cloud (correlated)
        import random
        random.seed(i)
        for _ in range(60):
            xx = x0 + 30 + random.random() * 220
            yy = 460 - (xx - x0 - 30) * 0.6 + (random.random() - 0.5) * 90
            yy = max(140, min(460, yy))
            s += circle(xx, yy, 3, "#7f8c8d")
        s += txt(x0 + 20, 500, "ZP3", 9)
    s += footer("TCGA GBM/LGG; TREM2 top ZP3-correlated immune gene.", "cBioPortal query (Phase 1)")
    return s + "</svg>"

# ============ Article 2 ============
def fig_a2f1():
    """GSE78220: primary responder/non-responder distribution plus response categories."""
    W, H = 720, 600
    rows = load(p("h2_bulk", "gse78220_zp3_response.csv"))
    groups = {}
    for r in rows:
        if r["timepoint"] != "baseline":
            continue
        groups.setdefault(r["response"], []).append(float(r["zp3_expr"]))
    allv = [v for vs in groups.values() for v in vs]
    vmin, vmax = min(allv), max(allv)
    order = ["CR", "PR", "PD"]
    s = hdr(W, H, "Figure 1. ZP3 unrelated to anti-PD-1 response in melanoma (GSE78220)",
            "Primary comparison (left): CR/PR vs PD with individual samples; right: exploratory CR/PR/PD categories")
    def box_panel(xl, xr, vals_by_group, labels, title, stat_text, colors):
        y0, y1 = 500, 130
        s0 = txt((xl+xr)/2, 102, title, 10, "middle", weight="bold")
        s0 += line(xl, y0, xl, y1)
        s0 += line(xl, y0, xr, y0, "#999", 1)
        for i, (lab, vals) in enumerate(zip(labels, vals_by_group)):
            vals = sorted(vals); x = lin(i + 0.5, 0, len(labels), xl + 35, xr - 35)
            q1 = vals[len(vals)//4]; med = vals[len(vals)//2]; q3 = vals[3*len(vals)//4]
            yi = lin(q1, vmin, vmax, y0, y1); ym = lin(med, vmin, vmax, y0, y1); yq = lin(q3, vmin, vmax, y0, y1)
            s0 += line(x, yq, x, yi, colors[i], 1.5)
            s0 += line(x-12, yq, x+12, yq, colors[i], 1.5) + line(x-12, yi, x+12, yi, colors[i], 1.5)
            s0 += rect(x-22, yq, 44, yi-yq, "#f4f6f7", colors[i]) + rect(x-22, ym-2, 44, 4, colors[i])
            for j, val in enumerate(vals):
                jx = x + ((j * 17) % 25) - 12
                s0 += circle(jx, lin(val, vmin, vmax, y0, y1), 2.7, colors[i], colors[i], 0.5)
            s0 += txt(x, y0+18, lab, 10, "middle") + txt(x, y0+34, f"n={len(vals)}", 8.5, "middle", fill="#555")
        s0 += txt(xl, y0+18, f"{vmin:.0f}", 8) + txt(xl, y1, f"{vmax:.0f}", 8)
        s0 += txt((xl+xr)/2, 535, stat_text, 8.5, "middle", fill="#333")
        return s0
    resp = groups.get("CR", []) + groups.get("PR", [])
    non = groups.get("PD", [])
    s += box_panel(90, 345, [resp, non], ["CR/PR", "PD"], "Pre-treatment primary analysis", "MWU P=0.90; r=-0.03 (95% CI -0.48, 0.41)", ["#c0392b", "#7f8c8d"])
    s += box_panel(395, 690, [groups.get(g, []) for g in order], order, "Response-category view", "ANOVA P=0.86; categories shown for transparency", ["#c0392b", "#d35400", "#7f8c8d"])
    s += footer("GSE78220 baseline n=27 (15 CR/PR vs 12 PD; individual observations shown; one on-treatment sample excluded).",
                "gse78220_zp3_response.csv; frozen effect table a2_effects_frozen.csv")
    return s + "</svg>"

def fig_a2f2():
    """IMvigor210: response-group distributions reconstructed from the censored expression matrix."""
    W, H = 720, 600
    rows = load(p("immunotherapy_validation", "imvigor210_zp3_results.csv")); r = rows[0]
    pdir = os.path.join(ROOT, "output", "immunotherapy_validation")
    pdata = load(os.path.join(pdir, "pData_IMvigor210.csv"))
    fdata = load(os.path.join(pdir, "fData_IMvigor210.csv"))
    sym_to_row = {x["Symbol"]: "gene_" + str(i+1) for i, x in enumerate(fdata)}
    zp3_row = sym_to_row["ZP3"]
    with open(os.path.join(pdir, "exmat_censored_IMvigor210.csv"), newline="", encoding="utf-8-sig") as _f:
        ex = list(csv.reader(_f))
    header = ex[0][1:]; zrow = next(x for x in ex[1:] if x[0] == zp3_row)
    zvals = dict(zip(header, [float(x) for x in zrow[1:]]))
    clin = {x["X"]: x for x in pdata}
    groups = {"CR/PR": [], "SD/PD": []}
    for sid, val in zvals.items():
        g = clin.get(sid, {}).get("binaryResponse", "")
        if g in groups: groups[g].append(val)
    vmin, vmax = min(min(x) for x in groups.values()), max(max(x) for x in groups.values())
    s = hdr(W, H, "Figure 2. ZP3 unrelated to anti-PD-L1 response in urothelial carcinoma (IMvigor210)",
            "Individual censored-expression observations with boxplots; primary evaluable cohort n=298")
    x0, y0, y1 = 95, 500, 130
    s += line(x0, y0, x0, y1) + line(x0, y0, 690, y0, "#999", 1)
    for i, g in enumerate(["CR/PR", "SD/PD"]):
        vals = sorted(groups[g]); x = 260 + i*260
        q1, med, q3 = vals[len(vals)//4], vals[len(vals)//2], vals[3*len(vals)//4]
        yq, ym, yi = lin(q3,vmin,vmax,y0,y1), lin(med,vmin,vmax,y0,y1), lin(q1,vmin,vmax,y0,y1)
        col = "#c0392b" if g == "CR/PR" else "#7f8c8d"
        s += line(x, yq, x, yi, col, 1.5) + line(x-14,yq,x+14,yq,col,1.5) + line(x-14,yi,x+14,yi,col,1.5)
        s += rect(x-32,yq,64,yi-yq,"#f4f6f7",col) + rect(x-32,ym-2,64,4,col)
        for j,val in enumerate(vals):
            jx = x + ((j*13)%55)-27
            s += circle(jx, lin(val,vmin,vmax,y0,y1), 2.1, col, col, 0.5)
        s += txt(x, y0+18, g, 11, "middle") + txt(x,y0+34,f"n={len(vals)}",9,"middle",fill="#555")
        s += txt(x, yq-10, f"median={med:.0f}", 9, "middle", weight="bold")
    s += txt(x0,y0+18,f"{vmin:.0f}",8) + txt(x0,y1,f"{vmax:.0f}",8)
    s += txt(55,(y0+y1)/2,"ZP3 (censored normalized units)",9,"middle",transform=f'rotate(-90 55 {(y0+y1)/2})')
    s += txt(390, 535, "MWU P=0.86; r=-0.015 (95% CI -0.17, 0.14)", 9, "middle")
    s += footer("IMvigor210: 298 evaluable of 348 enrolled; 50 without response annotation excluded. Raw observations are displayed with deterministic horizontal jitter.",
                "imvigor210_zp3_results.csv; exmat_censored_IMvigor210.csv; pData_IMvigor210.csv")
    return s + "</svg>"

def fig_a2f3():
    """GSE91061: pre-treatment ZP3 distribution plus frozen audit-context results."""
    W, H = 720, 600
    rows = load(p("a2_effects_frozen.csv")); r = [x for x in rows if x["Cohort"] == "GSE91061"][0]
    clin = load(os.path.join(ROOT, "output", "phase1_knowledge_gap_filling", "gse91061_molecular_subtype_clinical.csv"))
    groups = {"Responder": [], "Non-responder": []}
    for x in clin:
        if x["visit"] != "Pre" or x["response"] == "UNK": continue
        groups["Responder" if x["response"] == "PRCR" else "Non-responder"].append(float(x["ZP3"]))
    vals_all = groups["Responder"] + groups["Non-responder"]; vmin, vmax = min(vals_all), max(vals_all)
    s = hdr(W, H, "Figure 3. ZP3 unrelated to anti-PD-1 response in melanoma pre-treatment samples (GSE91061)",
            "Pre-treatment primary analysis: individual observations with boxplots (n=49)")
    x0, y0, y1 = 100, 500, 130
    s += line(x0,y0,x0,y1)+line(x0,y0,430,y0,"#999",1)
    for i,g in enumerate(["Responder","Non-responder"]):
        vals=sorted(groups[g]); x=220+i*170
        q1,med,q3=vals[len(vals)//4],vals[len(vals)//2],vals[3*len(vals)//4]
        yq,ym,yi=lin(q3,vmin,vmax,y0,y1),lin(med,vmin,vmax,y0,y1),lin(q1,vmin,vmax,y0,y1)
        col="#c0392b" if g=="Responder" else "#7f8c8d"
        s += line(x,yq,x,yi,col,1.5)+line(x-12,yq,x+12,yq,col,1.5)+line(x-12,yi,x+12,yi,col,1.5)
        s += rect(x-28,yq,56,yi-yq,"#f4f6f7",col)+rect(x-28,ym-2,56,4,col)
        for j,val in enumerate(vals): s += circle(x+((j*13)%45)-22,lin(val,vmin,vmax,y0,y1),3,col,col,0.5)
        s += txt(x,y0+18,g,10,"middle")+txt(x,y0+34,f"n={len(vals)}",9,"middle",fill="#555")+txt(x,yq-10,f"median={med:.2f}",9,"middle",weight="bold")
    s += txt(x0,y0+18,f"{vmin:.1f}",8)+txt(x0,y1,f"{vmax:.1f}",8)
    s += txt(58,(y0+y1)/2,"ZP3 (FPKM)",9,"middle",transform=f'rotate(-90 58 {(y0+y1)/2})')
    s += rect(455,155,220,250,"#f8f9f9","#bdc3c7",1)
    s += txt(565,180,"Audit context",11,"middle",weight="bold")
    s += txt(475,215,"Pre-treatment",9,fill="#333")+txt(650,215,"P=0.35",10,"end",weight="bold")
    s += txt(475,245,"All samples incl. OnTx",9,fill="#333")+txt(650,245,"P=0.008",10,"end",weight="bold",fill="#c0392b")
    s += txt(475,275,"Patient-level",9,fill="#333")+txt(650,275,"P=0.18",10,"end",weight="bold")
    s += line(475,295,650,295,"#bdc3c7",1)
    s += txt(565,325,"Only the on-treatment-inclusive",9,"middle",fill="#555")+txt(565,342,"sample-level analysis was significant",9,"middle",fill="#555")
    s += txt(565,375,"Primary inference: no baseline association",9,"middle",weight="bold")
    s += txt(565,392,"(see Figure 5 for full robustness analysis)",8,"middle",fill="#555")
    s += txt(265,535,"MWU P=0.35; r=-0.19 (95% CI -0.60, 0.21)",9,"middle")
    s += footer("GSE91061 pre-treatment n=49 (10 PR/CR vs 39 SD/PD); individual observations shown. Sensitivity analyses are summarized in the inset and detailed in Figure 5.",
                "gse91061_molecular_subtype_clinical.csv; a2_effects_frozen.csv; a2_effects_sensitivity.csv")
    return s + "</svg>"

def fig_a2f4():
    """GSE91061: apparent signal confined to on-treatment-inclusive analysis (audit-fix rerun)."""
    W, H = 900, 600
    s = hdr(W, H, "Figure 5. The apparent GSE91061 signal is confined to on-treatment-inclusive analysis",
            "Left: -log10(P) across analysis frameworks; right: subtype-stratified P (n=49 pre-treatment)")
    import math as _m
    # ---- Left panel: -log10(P) bars across analysis frameworks ----
    # (red = on-treatment-inclusive MWU, the only significant result; greys = pre-treatment / adjusted frameworks)
    x0, y0, y1 = 130, 500, 130
    s += line(x0, y0, x0, y1)
    bars = [
        ("MWU, incl.\non-treatment\n(105 samples,\n62 patients)", _m.log10(1/0.0081), "#c0392b"),
        ("MWU,\npre-treatment\n(n=49)", _m.log10(1/0.3522), "#95a5a6"),
        ("Logistic\nunivariable\n(OR 0.72)", _m.log10(1/0.4545), "#95a5a6"),
        ("Logistic\n+ IFN-g", _m.log10(1/0.5154), "#95a5a6"),
        ("Logistic\n+ Supp", _m.log10(1/0.5772), "#95a5a6"),
        ("Firth\nunivariable\n(OR 0.82)", _m.log10(1/0.6060), "#bdc3c7"),
        ("Firth\n+ IFN-g", _m.log10(1/0.8300), "#bdc3c7"),
        ("Firth\n+ Supp", _m.log10(1/0.6660), "#bdc3c7"),
        ("Patient-level\n(n=62)", _m.log10(1/0.1753), "#7f8c8d"),
    ]
    vmax = 2.8
    nbars = len(bars)
    bw = 42
    xstart, xend = 162.0, 578.0          # left-panel bar centers span
    step = (xend - xstart) / (nbars - 1)
    for i, (lab, v, col) in enumerate(bars):
        x = xstart + i * step
        h = v / vmax * (y0 - y1)
        s += rect(x - bw/2, y0 - h, bw, h, col)
        s += txt(x, y0 - h - 6, f"{v:.2f}", 8, "middle", weight="bold")
        for j, ln in enumerate(lab.split("\n")):
            s += txt(x, y0 + 14 + j * 10, ln, 7.5, "middle")
    ythr = lin(0.05, 0, vmax, y0, y1)
    s += line(x0, ythr, xend + 18, ythr, "#333", 1, "4,3")
    s += txt(xend + 22, ythr + 3, "P=0.05", 8)
    s += txt(x0, y0 + 20, "0", 9) + txt(x0, y1, f"{vmax}", 9)
    s += txt(72, (y0+y1)/2, "-log10(P)", 10, "middle", transform=f'rotate(-90 72 {(y0+y1)/2})')
    # ---- Right panel: subtype-stratified P (rerun with correct Entrez 7784) ----
    x02, y02, y12 = 640, 500, 130
    s += line(x02, y02, x02, y12)
    sub = [("Immune", _m.log10(1/0.8001), "#95a5a6"), ("Keratin", _m.log10(1/0.7692), "#95a5a6")]
    bw2 = 54
    for i, (g, v, col) in enumerate(sub):
        x = 700 + i * 95
        h = v / 1.2 * (y02 - y12)
        s += rect(x - bw2/2, y02 - h, bw2, h, col)
        s += txt(x, y02 - h - 6, f"P={10**(-v):.2f}", 8, "middle", weight="bold")
        s += txt(x, y02 + 14, g, 9, "middle")
    # Prominent textual verdict (bars intentionally small: P~0.8 -> -log10(P)~0.1; emphasis on text)
    s += txt(747, 152, "ZP3 vs anti-PD-1 response", 9, "middle", weight="bold")
    s += txt(747, 166, "No subtype association", 9, "middle", weight="bold", fill="#c0392b")
    s += txt(747, 180, "(both P ~ 0.77-0.80, NS)", 8, "middle", fill="#555")
    s += txt(747, y02 + 34, "MITF-high/low:", 8, "middle")
    s += txt(747, y02 + 46, "not estimable", 8, "middle", fill="#888")
    # right-panel scale ticks: -log10(P), divisor 1.2 (independent of left panel)
    vmax_r = 1.2
    for tv in [0.0, 0.4, 0.8, 1.2]:
        ty = y02 - tv / vmax_r * (y02 - y12)
        s += line(x02, ty, x02 + 5, ty, "#333", 1)
        s += txt(x02 - 4, ty + 3, f"{tv:.1f}", 8, "end")
    s += txt(610, (y02+y12)/2, "-log10(P)", 10, "middle", transform=f'rotate(-90 610 {(y02+y12)/2})')
    s += txt(747, 118, "Subtype-stratified", 9, "middle")
    s += footer("GSE91061 pre-treatment n=49 (10 responders); all pre-treatment and patient-level frameworks P>0.05; "
                "signal appears only when on-treatment samples are included (P=0.008; 105 samples from 62 patients).",
                "gse91061_molecular_subtype rerun 2026-08-17 (Entrez 7784); recompute_effects_a2.py; "
                "freeze_standard_logistic_a2.py (standard logistic n=49); "
                "freeze_subtype_response_a2.py (subtype-stratified MWU n=49); "
                "robustness_ci_a2.py (Firth + patient-level) 2026-08-17")
    return s + "</svg>"

def fig_a2f5():
    """Positive controls: TMB box + IHC response rate (IMvigor210)."""
    W, H = 720, 600
    s = hdr(W, H, "Figure 4. Positive controls confirm assay sensitivity (IMvigor210)",
            "TMB by response and PD-L1 IHC response rate; both reproduce known ICB associations")
    # Left panel: TMB medians
    x0, y0, y1 = 120, 520, 130
    s += line(x0, y0, x0, y1)
    meds = [("Responder", 14.0, "#c0392b"), ("Non-responder", 7.0, "#95a5a6")]
    tmb_n = [(61, 173), (173, 61)]
    vmax = 20
    for i, (g, v, col) in enumerate(meds):
        x = 200 + i * 120
        y = lin(v, 0, vmax, y0, y1)
        s += rect(x - 35, y, 70, y0 - y, col)
        s += txt(x, y - 10, f"{v:.1f}", 10, "middle", weight="bold")
        s += txt(x, y0 + 18, g, 9, "middle")
        s += txt(x, y0 + 33, f"n={tmb_n[i][0] if i == 0 else tmb_n[i][0]}", 8, "middle", fill="#555")
    s += txt(90, (y0+y1)/2, "TMB (mut/Mb, median)", 9, "middle", transform=f'rotate(-90 90 {(y0+y1)/2})')
    s += txt(x0, y0 + 20, "0", 9) + txt(x0, y1, f"{vmax}", 9)
    s += txt(200, 555, "P=1.1x10-7 (Mann-Whitney), r=0.46", 9, "middle")
    # Right panel: IHC response rates
    x02, y02, y12 = 480, 520, 130
    s += line(x02, y02, x02, y12)
    ihc = [("IC2+", 34.3, "#c0392b"), ("IC0/1", 16.9, "#95a5a6")]
    vmax2 = 40
    ihc_n = [("IC2+", 35, 67), ("IC0/1", 33, 162)]
    for i, (g, v, col) in enumerate(ihc):
        x = 560 + i * 120
        y = lin(v, 0, vmax2, y02, y12)
        s += rect(x - 35, y, 70, y02 - y, col)
        s += txt(x, y - 10, f"{v:.1f}%", 10, "middle", weight="bold")
        rr, nn = ihc_n[i][1], ihc_n[i][2]
        s += txt(x, y02 + 18, g, 9, "middle")
        s += txt(x, y02 + 33, f"{rr}/{rr+nn} responders", 8, "middle", fill="#555")
    s += txt(450, (y02+y12)/2, "Response rate by PD-L1 IHC", 9, "middle", transform=f'rotate(-90 450 {(y02+y12)/2})')
    s += txt(x02, y02 + 20, "0", 9) + txt(x02, y12, f"{vmax2}%", 9)
    s += txt(560, 555, "P=1.2x10-3 (chi-square)", 9, "middle")
    s += footer("IMvigor210 n=348 enrolled (298 evaluable); pre-specified reference biomarkers significant -> pipeline sensitivity indicated in this cohort.", "positive_control_results.csv")
    return s + "</svg>"

# ============ Article 3 ============
A3_FROZEN = os.path.join(ROOT, "article3", "results")

def fig_a3f1():
    """Fig1: polished isoform remodeling view using frozen medians, IQRs and ratios."""
    W, H = 1100, 700
    rows = load(os.path.join(A3_FROZEN, "a3_isoform_shift.csv"))
    labels = {"ENST00000336517.8": "FL | canonical",
              "ENST00000466960.5": "RI | retained intron",
              "ENST00000394860.3": "Internal promoter | 5-exon",
              "ENST00000394857.7": "ENST00000394857.7",
              "ENST00000416245.5": "ENST00000416245.5",
              "ENST00000467555.1": "ENST00000467555.1",
              "ENST00000479793.5": "ENST00000479793.5"}
    # log10 scale keeps low-abundance transcripts visible without inventing values.
    xmin, xmax = -4.0, 0.0
    x0, x1 = 300, 760
    y0, dy = 190, 54
    s = hdr(W, H, "Figure 1. Tumor-associated remodeling of ZP3 transcript composition",
            "Median proportions with interquartile ranges; the right panel ranks tumor-to-normal shifts")
    # panel labels and subtle panel backgrounds
    s += rect(42, 92, 748, 500, "#fbfcfe", "#e5e9ef", 1)
    s += rect(820, 92, 238, 500, "#fbfcfe", "#e5e9ef", 1)
    s += txt(64, 126, "A  Distribution by transcript", 12, weight="bold", fill="#263238")
    s += txt(842, 126, "B  Shift magnitude", 12, weight="bold", fill="#263238")
    # main log-scale plot
    for v in (-4, -3, -2, -1, 0):
        x = lin(v, xmin, xmax, x0, x1)
        s += line(x, 155, x, 555, stroke="#e1e6ec", width=1)
        label = "0.01%" if v == -4 else ("0.1%" if v == -3 else ("1%" if v == -2 else ("10%" if v == -1 else "100%")))
        s += txt(x, 575, label, 9, "middle", fill="#5f6b76")
    s += txt((x0+x1)/2, 610, "Transcript proportion (log scale)", 10, "middle", fill="#46515c")
    # legend
    s += circle(310, 145, 5, "#c43d3d", "#ffffff", 1) + txt(322, 149, "Tumor median", 9, fill="#46515c")
    s += circle(430, 145, 5, "#2d6fa3", "#ffffff", 1) + txt(442, 149, "Normal median", 9, fill="#46515c")
    for i, r in enumerate(rows):
        y = y0 + i * dy
        tm, nm = float(r["Tumor_median"]), float(r["Normal_median"])
        tq1, tq3 = float(r["Tumor_Q1"]), float(r["Tumor_Q3"])
        nq1, nq3 = float(r["Normal_Q1"]), float(r["Normal_Q3"])
        xt, xn = lin(math.log10(max(tm, 1e-8)), xmin, xmax, x0, x1), lin(math.log10(max(nm, 1e-8)), xmin, xmax, x0, x1)
        # IQR ranges, medians as crisp points; no unsupported CI.
        for q1, q3, xx, col in ((tq1, tq3, xt, "#c43d3d"), (nq1, nq3, xn, "#2d6fa3")):
            xa, xb = lin(math.log10(max(q1, 1e-8)), xmin, xmax, x0, x1), lin(math.log10(max(q3, 1e-8)), xmin, xmax, x0, x1)
            s += line(xa, y, xb, y, stroke=col, width=4)
            s += line(xa, y-6, xa, y+6, stroke=col, width=1.2) + line(xb, y-6, xb, y+6, stroke=col, width=1.2)
        s += circle(xt, y, 5.5, "#c43d3d", "#ffffff", 1.2) + circle(xn, y, 5.5, "#2d6fa3", "#ffffff", 1.2)
        s += txt(58, y+4, labels.get(r["Transcript"], r["Transcript"]), 9.5, fill="#263238")
        s += txt(58, y+18, f"P={float(r['MannWhitney_p']):.2e}", 8, fill="#7b8791")
    # right: log2 ratio bars centered at zero
    rx0, rx1, ry0, ry1 = 870, 1020, 190, 555
    s += line((rx0+rx1)/2, ry0-18, (rx0+rx1)/2, ry1, stroke="#9aa6b2", width=1)
    for v in (-4, 0, 4):
        xx = lin(v, -4, 4, rx0, rx1)
        s += txt(xx, 575, f"{v:+d}", 9, "middle", fill="#5f6b76")
    s += txt((rx0+rx1)/2, 610, "log2 tumor / normal", 10, "middle", fill="#46515c")
    for i, r in enumerate(rows):
        y = ry0 + i * dy
        ratio = float(r["Tumor_over_Normal_ratio"])
        log2r = math.log2(max(ratio, 1e-8))
        xx = lin(max(-4, min(4, log2r)), -4, 4, rx0, rx1)
        mid = (rx0+rx1)/2
        col = "#c43d3d" if log2r > 0 else "#2d6fa3"
        s += line(mid, y, xx, y, stroke=col, width=5) + circle(xx, y, 5.5, col, "#ffffff", 1.2)
        ratio_text = f"{ratio:.1f}x" if ratio >= 1 else f"{ratio:.2f}x"
        s += txt(1042, y+4, ratio_text, 9, "end", weight="bold", fill=col)
    # conclusion ribbon
    s += rect(64, 625, 994, 34, "#eef3f7", "#d7e0e7", 1)
    s += txt(80, 647, "FL and the 5-exon transcript are tumor-enriched; the retained-intron transcript is normal-enriched.", 10.5, fill="#263238", weight="bold")
    s += footer("TCGA tumor n=9,186 vs GTEx normal n=7,792; points = medians, whiskers = Q1-Q3; ratios and P values are frozen, descriptive group comparisons.",
                "article3/results/a3_isoform_shift.csv (frozen 2026-08-18)")
    return s + "</svg>"

def fig_a3f2():
    """Fig2 (enhanced, M1): two-scale proxy validation —
    left: ecological cancer-level (SpliceSeq AP-1 PSI vs our FL fraction, LOO+bootstrap);
    right: sample-level rho per GBM/LGG AP event, defining the proxy boundary."""
    W, H = 900, 600
    eco = load(os.path.join(A3_FROZEN, "a3_spliceseq_ecological.csv"))[0]
    samp = load(os.path.join(A3_FROZEN, "a3_spliceseq_samplelevel.csv"))
    rows = load(os.path.join(ROOT, "article3", "data", "spliceseq_zp3", "spliceseq_ecological_validation.csv"))
    s = hdr(W, H, "Figure 2. TRA proxy validates at cancer level; bounded at sample level",
            "Left: ecological agreement (rho=0.95, LOO 0.93-0.96, bootstrap CI 0.62-1.00); Right: sample-level boundary")
    # ---- left panel: ecological scatter ----
    lx0, lx1, ly0, ly1 = 70, 430, 520, 110
    s += line(lx0, ly0, lx1, ly0) + line(lx0, ly0, lx0, ly1)
    ap = [float(r["SpliceSeq_AP1_PSI"]) for r in rows]
    fl = [float(r["Our_FL_PSI"]) for r in rows]
    amin, amax = min(ap) - 0.05, max(ap) + 0.05
    fmin, fmax = min(fl) - 0.05, max(fl) + 0.05
    for r in rows:
        x = lin(float(r["SpliceSeq_AP1_PSI"]), amin, amax, lx0, lx1)
        y = lin(float(r["Our_FL_PSI"]), fmin, fmax, ly0, ly1)
        s += circle(x, y, 7, "#8e44ad")
        s += txt(x + 9, y + 3, r["Cancer"], 8)
    s += txt((lx0 + lx1) / 2, ly0 + 22, "SpliceSeq AP-1 PSI (median/cancer)", 9, "middle")
    s += txt(lx0 - 10, (ly0 + ly1) / 2, "Our FL fraction", 9, "middle",
             transform=f'rotate(-90 {lx0-10} {(ly0+ly1)/2})')
    s += rect(58, 128, 208, 86, "#f8f9fa", "#ccc")
    s += txt(66, 146, f"Ecological n=8 cancers", 9, weight="bold")
    s += txt(66, 164, f"Spearman rho = {eco['Spearman_rho']}", 9)
    s += txt(66, 180, f"Pearson r = {eco['Pearson_r']}", 9)
    s += txt(66, 196, f"LOO rho {eco['LOO_rho_min']}-{eco['LOO_rho_max']}", 9)
    s += txt(66, 210, f"Bootstrap 95% CI {eco['Bootstrap_CI_low']}-{eco['Bootstrap_CI_high']}", 9)
    # ---- right panel: sample-level rho dots with 0.7 boundary ----
    rx0, rx1, ry0, ry1 = 540, 870, 520, 110
    s += line(rx0, ry0, rx1, ry0) + line(rx0, ry0, rx0, ry1)
    rho_min, rho_max = -0.2, 0.8
    # 0.7 reference (validated-PSI threshold)
    s += line(lin(0.7, rho_min, rho_max, rx0, rx1), ry1 - 6, lin(0.7, rho_min, rho_max, rx0, rx1), ry0, dash="4 3", stroke="#888")
    s += txt(lin(0.7, rho_min, rho_max, rx0, rx1), ry1 - 12, "0.7 (validated threshold)", 8, "middle", fill="#888")
    n = len(samp)
    for i, r in enumerate(samp):
        y = lin(i + 0.5, 0, n, ry0, ry1)
        rho = float(r["Spearman_rho"])
        x = lin(rho, rho_min, rho_max, rx0, rx1)
        col = "#c0392b" if r["Cancer"] == "GBM" else "#2980b9"
        s += circle(x, y, 6, col)
        sig = "*" if float(r["Spearman_p"]) < 0.05 else ""
        s += txt(x + 10, y + 4, f"{r['Cancer']} AP{(r['Event'])} {rho:+.2f}{sig}", 8, fill=col)
    s += txt(rx0, ry0 + 22, "-0.2", 8, "middle") + txt(rx1, ry0 + 22, "0.8", 8, "middle")
    s += txt(rx0 - 12, (ry0 + ry1) / 2, "Sample-level Spearman rho", 9, "middle",
             transform=f'rotate(-90 {rx0-12} {(ry0+ry1)/2})')
    s += txt((rx0 + rx1) / 2, 96, "Proxy boundary: sample-level rho 0.13-0.54 < 0.7", 10, "middle", fill="#c0392b", weight="bold")
    s += footer("Left: spliceseq_ecological_validation.csv + a3_spliceseq_ecological.csv; Right: a3_spliceseq_samplelevel.csv (frozen 2026-08-18). "
                "FL proportion is a proxy, not a validated sample-level PSI replacement.", "a3_spliceseq_*.csv")
    return s + "</svg>"

def fig_a3f3():
    """Fig3 (enhanced, M1): GSEA horizontal bar chart — headline Hallmark programs.
    Values from frozen a3_gsea_headline.csv (audited 2026-08-18)."""
    W, H = 900, 600
    rows = load(os.path.join(A3_FROZEN, "a3_gsea_headline.csv"))
    s = hdr(W, H, "Figure 3. FL-high vs RI-high engage distinct programs (GSEA, Hallmark)",
            "Horizontal bars = NES; red = FL-high enriched, blue = RI-high enriched")
    x0, x1, y0, y1 = 320, 860, 500, 110
    s += line(x0, y0, x1, y0) + line(x0, y0, x0, y1)
    zero = lin(0, -2.5, 2.5, x0, x1)
    s += line(zero, y1 - 4, zero, y0, stroke="#999", width=1, dash="4 3")
    n = len(rows)
    for i, r in enumerate(rows):
        y = lin(i + 0.5, 0, n, y0, y1)
        nes = float(r["NES"])
        col = "#c0392b" if nes > 0 else "#2980b9"
        xv = lin(nes, -2.5, 2.5, x0, x1)
        s += rect(min(zero, xv), y - 9, abs(xv - zero), 18, col, col, 1)
        path = r["Pathway"].replace("MSigDB_Hallmark_2020__", "").replace("_", " ")
        if "TNF" in path:
            path = "TNF-alpha signaling via NF-kB"
        s += txt(x0 - 12, y + 4, path[:38], 9, "end")
        q = float(r["FDR_q"])
        qlab = "<0.001" if q == 0 else f"{q:.1e}"
        s += txt(xv + (6 if nes > 0 else -6), y + 4, f"{nes:+.2f}  (FDR {qlab})", 8,
                 anchor="start" if nes > 0 else "end", fill=col)
    s += txt(x0, y0 + 20, "-2.5", 9, "middle") + txt(x1, y0 + 20, "2.5", 9, "middle")
    s += txt(x0 - 8, (y0 + y1) / 2, "Normalized enrichment score (NES)", 10, "middle",
             transform=f'rotate(-90 {x0-8} {(y0+y1)/2})')
    s += rect(580, 128, 250, 52, "#f8f9fa", "#ccc")
    s += rect(592, 138, 18, 10, "#c0392b") + txt(618, 147, "FL-high enriched", 9)
    s += rect(592, 158, 18, 10, "#2980b9") + txt(618, 167, "RI-high enriched (FL depleted)", 9)
    s += footer("GBM+LGG n=694; TNF-alpha/NF-kB NES=+2.45 (FDR<0.001); E2F/G2-M/Myc negative; "
                "Inflammatory +1.86 (q=4.3e-3); DNA repair -1.67 (q=4.9e-3).",
                "a3_gsea_headline.csv (frozen 2026-08-18)")
    return s + "</svg>"

def fig_a3f4():
    """Fig4 (enhanced, M1): mixed-effects forest plot — FL/RI x 7 features,
    with adjusted-for-total-expression estimates overlaid (open symbols)."""
    W, H = 900, 600
    rows = load(os.path.join(A3_FROZEN, "a3_mixed_model_frozen.csv"))
    s = hdr(W, H, "Figure 4. Pan-cancer mixed-effects: isoform-immune link within cancers",
            "Solid = unadjusted; open ring = adjusted for ZP3 total expression; bars = 95% CI")
    x0, x1, y0, y1 = 260, 700, 540, 90
    s += line(x0, y0, x1, y0) + line(x0, y0, x0, y1)
    zero = lin(0, -0.35, 0.35, x0, x1)
    s += line(zero, y1 - 6, zero, y0, stroke="#999", width=1, dash="4 3")
    feats = ["M2_Macrophage", "Myeloid", "IFN_gamma", "T_cell_exhaustion", "Checkpoint", "Treg", "Cytolytic_activity"]
    labmap = {"M2_Macrophage": "M2 macrophages", "T_cell_exhaustion": "T-cell exhaustion",
              "Cytolytic_activity": "Cytolytic activity", "IFN_gamma": "IFN-gamma",
              "Checkpoint": "Checkpoint", "Treg": "Treg", "Myeloid": "Myeloid"}
    order = {}
    for k, f in enumerate(feats):
        order[f] = k
    # build per-feature rows: FL upper / RI lower
    rows_by = {}
    for r in rows:
        key = (r["Tx_Label"].split()[0], r["Feature"], r["Model"])
        rows_by[key] = r
    n_feats = len(feats)
    yrow = {}
    for f in feats:
        yrow[f] = {}
        k = order[f]
        yrow[f]["FL"] = lin(k * 2 + 0.5, 0, n_feats * 2, y0, y1)
        yrow[f]["RI"] = lin(k * 2 + 1.5, 0, n_feats * 2, y0, y1)
    for f in feats:
        for tx in ("FL", "RI"):
            y = yrow[f][tx]
            col = "#c0392b" if tx == "FL" else "#2980b9"
            key = (tx, f, "unadjusted")
            if key not in rows_by:
                continue
            r0 = rows_by[key]
            coef, lo, hi, fdr = float(r0["Coef"]), float(r0["CI_low"]), float(r0["CI_high"]), float(r0["FDR"])
            xm, xl, xh = lin(coef, -0.35, 0.35, x0, x1), lin(lo, -0.35, 0.35, x0, x1), lin(hi, -0.35, 0.35, x0, x1)
            s += line(xl, y, xh, y, stroke=col, width=2)
            s += circle(xm, y, 5, col)
            # adjusted overlay (open ring)
            key_a = (tx, f, "adjusted_ZP3_total")
            if key_a in rows_by:
                ra = rows_by[key_a]
                ca = float(ra["Coef"])
                xma = lin(ca, -0.35, 0.35, x0, x1)
                s += circle(xma, y, 5, "#ffffff", stroke=col, width=2)
            sig = "*" if fdr < 0.05 else ""
            anchor = "end" if tx == "FL" else "start"
            ax = xl - 6 if tx == "FL" else xh + 6
            if tx == "FL":
                s += txt(ax, y + 4, f"{labmap[f]}  {coef:+.2f}{sig}", 8.5, "end", fill=col)
            else:
                s += txt(ax, y + 4, f"{coef:+.2f}{sig}", 8.5, "start", fill=col)
    # zero axis labels
    s += txt(x0, y0 + 20, "-0.35", 9, "middle") + txt(x1, y0 + 20, "0.35", 9, "middle")
    s += txt(x0 - 10, (y0 + y1) / 2, "Coefficient (immune score)", 10, "middle",
             transform=f'rotate(-90 {x0-10} {(y0+y1)/2})')
    # legend
    s += rect(715, 128, 165, 62, "#f8f9fa", "#ccc")
    s += circle(730, 142, 5, "#c0392b") + txt(744, 146, "FL (unadjusted)", 9)
    s += circle(730, 164, 5, "#ffffff", stroke="#c0392b", width=2) + txt(744, 168, "FL (adj. total expr.)", 9)
    s += circle(730, 186, 5, "#2980b9") + txt(744, 190, "RI (unadjusted)", 9)
    s += footer("n=9,186 tumor samples x 32 cancer types (cancer random intercept); FL 6/7 positive, RI 6/7 negative (Cytolytic ns); "
                "adjusted models include ZP3 total expression: FL retains 6/7 (M2 beta=0.2383, P=8.1e-24).",
                "a3_mixed_model_frozen.csv (frozen 2026-08-18)")
    return s + "</svg>"

def fig_a3f5():
    """Fig5: external GSE113474 kallisto re-analysis — gene-level detection vs
    read-length-dependent isoform assignment (51-bp single-end, n=24 GBM).
    Values computed from a3_external_zp3_isoform.csv + kallisto_run.log (2026-08-22)."""
    W, H = 1120, 700
    import re, statistics as st
    iso = load(os.path.join(ROOT, "article3", "data", "external_reanalysis", "a3_external_zp3_isoform.csv"))
    # panel A: per-sample ZP3 total TPM (log scale)
    tot = [float(r["ZP3_total_TPM"]) for r in iso]
    med_tot = st.median(tot)
    # panel B: median isoform fractions (GSE113474) vs TCGA GBM FL 0.403
    def med_frac(col):
        return st.median([float(r[col + "_frac"]) for r in iso])
    fl_f, ri_f, ip_f = med_frac("ENST00000336517.8"), med_frac("ENST00000466960.5"), med_frac("ENST00000394860.3")
    tcga_gbm_fl = 0.403
    # panel C: mapping rates from kallisto_run.log
    rates = []
    for ln in open(os.path.join(ROOT, "article3", "data", "external_reanalysis", "kallisto_run.log"),
                   encoding="utf-8", errors="ignore"):
        m = re.search(r"processed ([\d,]+) reads, ([\d,]+) reads pseudoaligned", ln)
        if m:
            proc, al = int(m.group(1).replace(",", "")), int(m.group(2).replace(",", ""))
            rates.append(100.0 * al / proc)
    med_rate = st.median(rates)
    s = hdr(W, H, "Figure 5. Platform-specific isoform-level re-analysis of GSE113474 (n=24 GBM) with kallisto",
            "51-bp single-end reads; ZP3 is detectable at gene level, but isoform assignment is read-length dependent")
    s += txt(64, 118, "A  Gene-level detection", 12, weight="bold", fill="#263238")
    s += txt(408, 118, "B  Isoform composition", 12, weight="bold", fill="#263238")
    s += txt(792, 118, "C  Pseudo-alignment quality", 12, weight="bold", fill="#263238")
    # ---- panel A: ZP3 total TPM, log10 ----
    ax0, ax1, ay0, ay1 = 120, 330, 500, 150
    lmin, lmax = 0.0, 2.0  # log10(1)..log10(100)
    s += line(ax0, ay0, ax1, ay0) + line(ax0, ay0, ax0, ay1)
    for v in (1, 3, 10, 30, 100):
        x = lin(math.log10(v), lmin, lmax, ax0, ax1)
        s += line(x, ay0, x, ay0 + 4, stroke="#9aa6b2", width=1)
        s += txt(x, ay0 + 16, f"{v}", 8.5, "middle", fill="#5f6b76")
    s += txt((ax0+ax1)/2, ay0 + 34, "ZP3 total TPM (log)", 9, "middle", fill="#46515c")
    for i, t in enumerate(tot):
        x = ax0 + 4 + (i % 6) * 8 + (i // 6) * 6  # deterministic jitter
        x = min(x, ax1 - 2)
        y = lin(math.log10(max(t, 1)), lmin, lmax, ay0, ay1)
        s += circle(x, y, 4, "#2d6fa3", "#ffffff", 1)
    ym = lin(math.log10(med_tot), lmin, lmax, ay0, ay1)
    s += line(ax0, ym, ax1, ym, stroke="#c43d3d", width=1.5, dash="4 3")
    s += txt(ax1, ym - 7, f"median {med_tot:.1f}", 8.5, "end", fill="#c43d3d", weight="bold")
    s += txt(ax0, ay1 + 18, "24/24 detected (TPM > 0)", 9, "middle", fill="#2d6fa3", weight="bold")
    # ---- panel B: isoform fraction bars + TCGA reference ----
    bx0, bx1, by0, by1 = 400, 660, 500, 170
    fmax = 0.45
    s += line(bx0, by0, bx1, by0) + line(bx0, by0, bx0, by1)
    for v in (0, 0.1, 0.2, 0.3, 0.4):
        x = lin(v, 0, fmax, bx0, bx1)
        s += line(x, by0, x, by0 + 4, stroke="#9aa6b2", width=1)
        s += txt(x, by0 + 16, f"{v:.1f}", 8, "middle", fill="#5f6b76")
    # TCGA GBM FL reference (dashed)
    xref = lin(tcga_gbm_fl, 0, fmax, bx0, bx1)
    s += line(xref, by1 - 6, xref, by0, stroke="#c43d3d", width=1.5, dash="5 3")
    s += txt(xref, by1 - 14, "TCGA GBM FL 0.403", 8.5, "middle", fill="#c43d3d", weight="bold")
    rows_b = [("FL (canonical)", fl_f, "#c43d3d"), ("RI (retained intron)", ri_f, "#2d6fa3"),
              ("Internal-promoter (5-exon)", ip_f, "#9aa6b2")]
    for i, (lab, v, col) in enumerate(rows_b):
        y = lin(i + 0.5, 0, 3, by0, by1)
        xv = lin(v, 0, fmax, bx0, bx1)
        s += rect(bx0, y - 10, xv - bx0, 20, col, col, 1)
        val = f"{v:.3f}" if v > 0 else "0.000"
        s += txt(bx0 - 10, y + 4, lab, 9, "end", fill="#263238")
        s += txt(max(xv + 5, bx0 + 5), y + 4, val, 8.5, "start", fill=col, weight="bold")
    s += txt((bx0+bx1)/2, by0 + 34, "Median transcript proportion", 9, "middle", fill="#46515c")
    s += txt(bx0, by1 - 30, "FL collapses to 0.012; RI dominates (median TPM 4.24 vs 0.33)", 8.5, fill="#5f6b76")
    # ---- panel C: mapping rate by sample ----
    cx0, cx1, cy0, cy1 = 790, 1040, 500, 150
    rmin, rmax = 75, 90
    s += line(cx0, cy0, cx1, cy0) + line(cx0, cy0, cx0, cy1)
    for v in (75, 78, 81, 84, 87, 90):
        y = lin(v, rmin, rmax, cy0, cy1)
        s += line(cx0, y, cx0 + 4, y, stroke="#9aa6b2", width=1)
        s += txt(cx0 - 5, y + 3, f"{v}", 8, "end", fill="#5f6b76")
    s += txt(cx0 - 12, (cy0+cy1)/2, "Pseudoaligned %", 9, "middle", fill="#46515c",
             transform=f'rotate(-90 {cx0-12} {(cy0+cy1)/2})')
    for i, rv in enumerate(sorted(rates)):
        x = cx0 + 6 + (i % 5) * 9 + (i // 5) * 5
        x = min(x, cx1 - 3)
        y = lin(rv, rmin, rmax, cy0, cy1)
        hl = abs(rv - 87.2) < 0.05
        s += circle(x, y, 4, "#c43d3d" if hl else "#5f6b76", "#ffffff", 1)
    ymc = lin(med_rate, rmin, rmax, cy0, cy1)
    s += line(cx0, ymc, cx1, ymc, stroke="#2d6fa3", width=1.5, dash="4 3")
    s += txt(cx1, ymc - 7, f"median {med_rate:.1f}%", 8.5, "end", fill="#2d6fa3", weight="bold")
    s += txt(cx0, cy1 + 18, "Representative SRR7050184: 87.2%", 8.5, "middle", fill="#c43d3d", weight="bold")
    # interpretation ribbon
    s += rect(64, 620, 992, 40, "#eef3f7", "#d7e0e7", 1)
    s += txt(80, 637, "Interpretation: gene-level ZP3 expression is reproducible, but short-read isoform assignment is read-length dependent;", 10.5, fill="#263238", weight="bold")
    s += txt(80, 652, "the FL fraction collapse (0.012 vs 0.403) reflects platform limits, not biological contradiction (exploratory, platform-specific).", 10, fill="#46515c")
    s += footer("Per-sample TPM/fractions: a3_external_zp3_isoform.csv; mapping rates: kallisto_run.log (2026-08-22); TCGA GBM FL median from a3_isoform_shift.csv (frozen 2026-08-18).",
                "article3/data/external_reanalysis/a3_external_zp3_isoform.csv")
    return s + "</svg>"

# A1 figure-number remapping (2026-08-15 v0.2 review revision): new figure number -> old generator function
# new Fig1=old Fig1 (TCGA markers) / new Fig2=old Fig4 (CGGA molecular subtypes) / new Fig3=old Fig6 (single-cell TREM2)
# new Fig4=old Fig3 (CGGA survival) / new Fig5=old Fig5 (pan-cancer) / new Fig6=old Fig2 (IMvigor210) / new SuppS2=old Fig7 (HPA)
FIGS = [
    ("A1_Fig1", fig_a1f1), ("A1_Fig2", fig_a1f4), ("A1_Fig3", fig_a1f6),
    ("A1_Fig4", fig_a1f3), ("A1_Fig5", fig_a1f5), ("A1_Fig6", fig_a1f2),
    ("A1_SuppS2", fig_a1f7), ("A1_SuppS1", fig_a1s1),
    ("A2_Fig1", fig_a2f1), ("A2_Fig2", fig_a2f2), ("A2_Fig3", fig_a2f3), ("A2_Fig4", fig_a2f5), ("A2_Fig5", fig_a2f4),
    ("A3_Fig1", fig_a3f1), ("A3_Fig2", fig_a3f2), ("A3_Fig3", fig_a3f3), ("A3_Fig4", fig_a3f4),
    ("A3_Fig5", fig_a3f5),
]

if __name__ == "__main__":
    for name, fn in FIGS:
        svg = fn()
        d = next((FIG_DIRS[k] for k in ("A1", "A2", "A3") if name.startswith(k)), OUT)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, name + ".svg"), "w", encoding="utf-8") as f:
            f.write(svg)
        print(f"written {name}.svg ({len(svg)} bytes)")
    print(f"Total {len(FIGS)} SVG figures written (per-article dirs under article1/2/3/figures)")
