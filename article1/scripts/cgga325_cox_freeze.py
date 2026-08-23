# -*- coding: utf-8 -*-
"""
CGGA-325 Cox 生存分析结果固化（2026-08-14 全流程审计要求）。
背景：cgga325_validation.py 第 162-168 行 Cox 仅 stdout 打印、无落盘；
本次重跑并将 HR/CI/p 写入 CSV，供图注骨架 A1 Fig3 引用。
"""
import os, pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from scipy import stats

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output", "cgga_validation")
CLIN = os.path.join(BASE, "CGGA.mRNAseq_325_clinical.20200506.txt")
EXPR = os.path.join(BASE, "CGGA.mRNAseq_325.RSEM-genes.20200506.txt")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article1", "results", "cgga325_cox_results.csv")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

clin = pd.read_csv(CLIN, sep="\t")
expr = pd.read_csv(EXPR, sep="\t", index_col=0)
zp3 = expr.loc["ZP3"]
common = sorted(set(clin["CGGA_ID"]) & set(zp3.index))
ca = clin[clin["CGGA_ID"].isin(common)].set_index("CGGA_ID").loc[common]
za = zp3[common]

# --- 生存分析 ---
surv = pd.DataFrame({
    "time": ca["OS"].values,
    "event": ca["Censor (alive=0; dead=1)"].values,
    "zp3": za.values,
}, index=ca.index).dropna()

# Log-rank（中位分组）
med = surv["zp3"].median()
high = surv[surv["zp3"] >= med]
low = surv[surv["zp3"] < med]
lr = logrank_test(high["time"], low["time"],
                  event_observed_A=high["event"], event_observed_B=low["event"])

# Cox（连续 ZP3，单变量）
cph = CoxPHFitter()
cph.fit(surv[["time", "event", "zp3"]], duration_col="time", event_col="event")
hr = cph.hazard_ratios_["zp3"]
p = cph.summary["p"]["zp3"]
# 注意：lifelines confidence_intervals_ 返回 log-hazard (beta) 尺度 CI，需 exp 转换为 HR 尺度
import numpy as np
se_zp3 = cph.standard_errors_["zp3"]
beta_zp3 = cph.params_["zp3"]
ci_lo, ci_hi = float(np.exp(beta_zp3 - 1.96 * se_zp3)), float(np.exp(beta_zp3 + 1.96 * se_zp3))

# Cox（多变量：+年龄 +分级 +IDH +1p19q，若字段可用）
row_uni = {
    "Cohort": "CGGA-325", "Model": "Cox univariate (continuous ZP3)",
    "n": len(surv), "Events": int(surv["event"].sum()),
    "HR": hr, "CI_low": ci_lo, "CI_high": ci_hi, "P": p,
    "Logrank_p_median_split": lr.p_value,
    "ZP3_median_split": med,
    "High_n": len(high), "Low_n": len(low),
}

# 多变量：尝试调整 age/grade/IDH/1p19q
rows = [row_uni]
try:
    cols_avail = []
    colmap = {"Age": "Age", "Grade": "Grade", "IDH_mutation_status": "IDH",
              "1p19q_codeletion_status": "1p19q"}
    m_surv = surv.copy()
    for src, dst in colmap.items():
        if src in ca.columns:
            vals = ca.loc[surv.index, src].values
            if dst == "Grade":  # WHO II/III/IV -> 2/3/4
                vals = [{"WHO II": 2, "WHO III": 3, "WHO IV": 4}.get(str(v), None) for v in vals]
            elif dst == "IDH":  # Wildtype/Mutant -> 0/1
                vals = [{"Wildtype": 0, "Mutant": 1}.get(str(v), None) for v in vals]
            elif dst == "1p19q":  # Non-codel/Codel -> 0/1
                vals = [{"Non-codel": 0, "Codel": 1}.get(str(v), None) for v in vals]
            m_surv[dst] = vals
            cols_avail.append(dst)
    m_surv = m_surv.dropna()
    print(f"多变量 Cox: 可用样本 n={len(m_surv)} (协变量 {cols_avail})")
    if len(cols_avail) >= 1 and len(m_surv) >= 50:
        cph_m = CoxPHFitter()
        cph_m.fit(m_surv[["time", "event", "zp3"] + cols_avail],
                  duration_col="time", event_col="event")
        row_multi = {
            "Cohort": "CGGA-325", "Model": f"Cox multivariable (ZP3 + {','.join(cols_avail)})",
            "n": len(m_surv), "Events": int(m_surv["event"].sum()),
            "HR": cph_m.hazard_ratios_["zp3"],
            "CI_low": float(np.exp(cph_m.params_["zp3"] - 1.96 * cph_m.standard_errors_["zp3"])),
            "CI_high": float(np.exp(cph_m.params_["zp3"] + 1.96 * cph_m.standard_errors_["zp3"])),
            "P": cph_m.summary["p"]["zp3"],
            "Logrank_p_median_split": "", "ZP3_median_split": "", "High_n": "", "Low_n": "",
        }
        rows.append(row_multi)
    else:
        print("多变量 Cox 跳过：无可用协变量或样本不足")
except Exception as e:
    print("多变量 Cox 未完成:", e)

res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print("已固化:", OUT)
print(res.to_string())
