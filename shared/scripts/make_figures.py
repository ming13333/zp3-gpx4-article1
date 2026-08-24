# -*- coding: utf-8 -*-
"""H2/H3 figures: KM survival curves + immune marker correlation bar chart + ZP3-TREM2 scatter plot."""
import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"figure.dpi": 130, "font.size": 10, "axes.grid": True,
                     "grid.alpha": 0.3, "axes.axisbelow": True})

def load_patient_csv(name):
    df = pd.read_csv(os.path.join(BASE, name))
    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "patientId"}).set_index("patientId")
    return df

def get_p_from_log(study):
    """2026-08-10 recalculation fix: no longer parse from old logs (old logs had buggy 0.902/0.954),
    instead independently recompute log-rank p with standard hypergeometric variance (GBM≈0.353, LGG≈0.384)."""
    from scipy import stats as _st
    osf = "h2_gbm_tcga_zp3_os.csv" if "gbm" in study else "h2_lgg_tcga_zp3_os.csv"
    dd = pd.read_csv(os.path.join(BASE, osf)).dropna(subset=["ZP3", "time", "event"])
    d = dd["time"].values; e = dd["event"].values; g = (dd["ZP3"] > dd["ZP3"].median()).astype(int).values
    o = np.argsort(d); d, e, g = d[o], e[o], g[o]
    Oe = 0.0; V = 0.0; n1 = int((g == 1).sum()); n0 = int((g == 0).sum())
    for t in np.unique(d):
        at = (d == t)
        n1t = int(((g == 1) & at).sum()); n0t = int(((g == 0) & at).sum())
        dj = int(((g == 1) & at & (e == 1)).sum()) + int(((g == 0) & at & (e == 1)).sum())
        nt = n1 + n0
        if nt > 1 and dj > 0:
            Oe += int(((g == 1) & at & (e == 1)).sum()) - dj * n1 / nt
            V += n1 * n0 * dj * (nt - dj) / (nt * nt * (nt - 1))
        n1 -= n1t; n0 -= n0t
    chi2 = Oe * Oe / V if V > 0 else 0.0
    return f"{_st.chi2.sf(chi2, 1):.3f}"

# ---------------- H2: Kaplan-Meier (High vs Low ZP3) ----------------
def km_step(time, event, group):
    out = {}
    for g in [0, 1]:
        t = time[group == g]; e = event[group == g]
        if len(t) == 0:
            continue
        order = np.argsort(t); t, e = t[order], e[order]
        times = np.unique(t); S = 1.0; xs, ys = [0.0], [1.0]
        for ti in times:
            at_risk = np.sum(t >= ti); d = np.sum((t == ti) & (e == 1))
            if at_risk > 0 and d > 0:
                S *= (1 - d / at_risk)
            xs.append(ti); ys.append(S)
        out[g] = (xs, ys, len(t))
    return out

for study, csv in [("gbm_tcga", "h2_gbm_tcga_zp3_os.csv"),
                   ("lgg_tcga", "h2_lgg_tcga_zp3_os.csv")]:
    df = load_patient_csv(csv)
    t, ev, g = df["time"].values, df["event"].values, df["group"].values
    km = km_step(t, ev, g)
    p = get_p_from_log(study)
    plt.figure(figsize=(4.6, 3.8))
    colors = {1: "#c0392b", 0: "#2c7fb8"}
    labels = {1: "High ZP3", 0: "Low ZP3"}
    for grp in [0, 1]:
        if grp in km:
            xs, ys, n = km[grp]
            plt.step(xs, ys, where="post", color=colors[grp], lw=2,
                     label=f"{labels[grp]} (n={n})")
    plt.ylim(0, 1.02); plt.xlim(0, max(t) * 1.02)
    plt.xlabel("Overall survival (months)"); plt.ylabel("Survival probability")
    plt.title(f"{study}: OS by ZP3 (median split)\nlog-rank p={p}")
    plt.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(BASE, f"fig_km_{study}.png")); plt.close()
    print(f"KM {study} done (p={p})")

# ---------------- H3: correlation bar (GBM & LGG) ----------------
h3_gbm = pd.read_csv(os.path.join(BASE, "h3_gbm_tcga_zp3_immuno.csv"))
h3_lgg = pd.read_csv(os.path.join(BASE, "h3_lgg_tcga_zp3_immuno.csv"))
rg = h3_gbm.set_index("gene")["spearman_rho"]; pg = h3_gbm.set_index("gene")["p"]
rl = h3_lgg.set_index("gene")["spearman_rho"]; pl = h3_lgg.set_index("gene")["p"]
order = rg.abs().sort_values(ascending=False).index.tolist()

fig, axes = plt.subplots(1, 2, figsize=(9.2, 6.2))
for ax, rr, pp, title in [(axes[0], rg, pg, "GBM (n=158)"),
                          (axes[1], rl, pl, "LGG (n=512)")]:
    vals = [rr[g] for g in order]
    sig = [pp[g] < 0.05 for g in order]
    cols = ["#c0392b" if s else "#7f8c8d" for s in sig]
    y = np.arange(len(order))[::-1]
    ax.barh(y, vals, color=cols)
    ax.set_yticks(y); ax.set_yticklabels(order, fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Spearman rho (ZP3 vs gene)")
    ax.set_title(f"H3 {title}\nZP3 vs immune markers (red = nominal sig.)")
    if "TREM2" in order:
        yi = y[order.index("TREM2")]
        ax.text(rr["TREM2"] + (0.005 if rr["TREM2"] >= 0 else -0.005), yi,
                " TREM2", va="center", fontsize=8, fontweight="bold", color="#c0392b")
    ax.set_xlim(-0.2, 0.28)
fig.suptitle("H3: correlation of bulk ZP3 with immunosuppressive markers", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(os.path.join(BASE, "fig_h3_immuno_bar.png")); plt.close()
print("H3 bar done")

# ---------------- ZP3 vs TREM2 scatter ----------------
for study, exprcsv, oscsv, outname in [("GBM", "expr_gbm_tcga_patient.csv",
                                         "h2_gbm_tcga_zp3_os.csv", "fig_scatter_gbm"),
                                        ("LGG", "expr_lgg_tcga_patient.csv",
                                         "h2_lgg_tcga_zp3_os.csv", "fig_scatter_lgg")]:
    expr = load_patient_csv(exprcsv); osdf = load_patient_csv(oscsv)
    common = expr.index.intersection(osdf.index)
    zp3 = osdf.loc[common, "ZP3"].values
    if "TREM2" not in expr.columns:
        continue
    trem = expr.loc[common, "TREM2"].values
    m = ~np.isnan(zp3) & ~np.isnan(trem)
    zp3, trem = zp3[m], trem[m]
    r, p = stats.spearmanr(zp3, trem)
    plt.figure(figsize=(4.4, 4.0))
    plt.scatter(zp3, trem, s=12, alpha=0.5, color="#2c7fb8")
    plt.xlabel("ZP3 (RSEM)"); plt.ylabel("TREM2 (RSEM)")
    plt.title(f"{study}: ZP3 vs TREM2\nrho={r:.3f}, p={p:.2e}, n={len(zp3)}")
    plt.tight_layout()
    plt.savefig(os.path.join(BASE, outname + ".png")); plt.close()
    print(f"scatter {study} done r={r:.3f} p={p:.2e}")
print("ALL FIGURES DONE")
