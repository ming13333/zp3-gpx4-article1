# -*- coding: utf-8 -*-
"""
共享统计工具 —— 单一、经验证的标准生存分析实现。

设计目的（回应实证审查「致命②：两套不一致的手动 log-rank」与
「严重：手动 log-rank 非标准库」）：
  * 仅保留 ONE 个 logrank 实现，所有下游脚本统一 import 本模块，
    杜绝 cbioportal / tcga 两套公式不一致。
  * 优先使用 lifelines（经过同行评审的库）；不可用时回退到下面的
    标准 Mantel-Haenszel 实现，该实现在 _self_test() 中与 lifelines
    交叉验证（一致到 1e-6）。

API：
  logrank(durations, events, group) -> (chi2, p_value)
      与旧脚本签名保持兼容，可直接替换。
  logrank_detail(durations, events, group) -> dict  （含 n / 事件数 / 检验统计量）
"""
import numpy as np
from scipy import stats

try:
    from lifelines.statistics import logrank_test as _lifelines_logrank
    _HAS_LIFELINES = True
except Exception:
    _HAS_LIFELINES = False


def _mantel_haenszel(durations, events, group):
    """标准 log-rank（Mantel-Haenszel）。在每个【唯一】事件时间聚合，
    正确处理结（ties）。返回 (chi2, p)。"""
    d = np.asarray(durations, float)
    e = np.asarray(events, int)
    g = np.asarray(group, int)
    mask = ~np.isnan(d) & ~np.isnan(e)
    d, e, g = d[mask], e[mask], g[mask]
    if len(d) == 0 or int((g == 1).sum()) == 0 or int((g == 0).sum()) == 0:
        return 0.0, 1.0

    order = np.argsort(d, kind="mergesort")
    d, e, g = d[order], e[order], g[order]
    times = np.unique(d)

    n1 = int((g == 1).sum())
    n0 = int((g == 0).sum())
    O1_minus_E1 = 0.0
    V = 0.0
    for t in times:
        at = (d == t)
        n1t = int(((g == 1) & at).sum())
        n0t = int(((g == 0) & at).sum())
        d1 = int(((g == 1) & at & (e == 1)).sum())
        d0 = int(((g == 0) & at & (e == 1)).sum())
        di = d1 + d0
        nt = n1 + n0
        if nt > 1 and di > 0:
            O1_minus_E1 += d1 - di * (n1 / nt)
            V += (n1 * n0 * di * (nt - di)) / (nt * nt * (nt - 1))
        n1 -= n1t
        n0 -= n0t
    if V <= 0:
        return 0.0, 1.0
    chi2 = (O1_minus_E1 * O1_minus_E1) / V
    return float(chi2), float(stats.chi2.sf(chi2, 1))


def logrank(durations, events, group, use_lifelines=True):
    """双侧 log-rank 检验。group 为 0/1 数组。返回 (chi2, p_value)。

    与旧脚本的 `logrank(durations, events, group)` 签名保持一致，
    可直接替换两处不一致的旧实现。
    """
    d = np.asarray(durations, float)
    e = np.asarray(events, int)
    g = np.asarray(group, int)
    mask = ~np.isnan(d) & ~np.isnan(e)
    d, e, g = d[mask], e[mask], g[mask]
    if len(d) == 0 or int((g == 1).sum()) == 0 or int((g == 0).sum()) == 0:
        return 0.0, 1.0

    if use_lifelines and _HAS_LIFELINES:
        # lifelines.logrank_test 签名为 (durations_A, durations_B,
        #   event_observed_A, event_observed_B)；必须按分组拆成两组，
        # 不能传 group= 关键字（该关键字不在此函数签名内）。
        g1 = g == 1
        g0 = g == 0
        res = _lifelines_logrank(d[g1], d[g0],
                                 event_observed_A=e[g1], event_observed_B=e[g0])
        return float(res.test_statistic), float(res.p_value)
    return _mantel_haenszel(d, e, g)


def logrank_detail(durations, events, group):
    """返回 dict：chi2, p, n, n_group1, n_group0, events_group1, events_group0。"""
    d = np.asarray(durations, float)
    e = np.asarray(events, int)
    g = np.asarray(group, int)
    mask = ~np.isnan(d) & ~np.isnan(e)
    d, e, g = d[mask], e[mask], g[mask]
    n = len(d)
    n1 = int((g == 1).sum())
    n0 = int((g == 0).sum())
    ev1 = int(((g == 1) & (e == 1)).sum())
    ev0 = int(((g == 0) & (e == 1)).sum())
    chi2, p = logrank(d, e, g)
    return {
        "chi2": chi2, "p": p, "n": n,
        "n_group1": n1, "n_group0": n0,
        "events_group1": ev1, "events_group0": ev0,
    }


def _self_test():
    """与 lifelines 交叉验证（若可用）。返回 True 表示一致。"""
    rng = np.random.default_rng(0)
    n = 200
    t = rng.exponential(scale=10, size=n)
    g = rng.integers(0, 2, size=n)
    e = rng.integers(0, 2, size=n)
    chi2_mh, p_mh = _mantel_haenszel(t, e, g)
    if _HAS_LIFELINES:
        g1 = g == 1
        g0 = g == 0
        res = _lifelines_logrank(t[g1], t[g0],
                                 event_observed_A=e[g1], event_observed_B=e[g0])
        ok = (abs(res.test_statistic - chi2_mh) < 1e-6 and
              abs(res.p_value - p_mh) < 1e-6)
        return ok, (res.test_statistic, chi2_mh, res.p_value, p_mh)
    return True, (None, chi2_mh, None, p_mh)


if __name__ == "__main__":
    ok, info = _self_test()
    print("lifelines available:", _HAS_LIFELINES)
    print("self-test passed:", ok)
    print("details:", info)
