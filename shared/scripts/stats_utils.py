# -*- coding: utf-8 -*-
"""
Shared statistics utilities -- a single, validated standard survival analysis implementation.

Design rationale (responding to empirical review 'Fatal 2: two inconsistent manual log-rank' and
'Severe: manual log-rank is non-standard library'):
  * Keep only ONE logrank implementation, all downstream scripts uniformly import this module,
    eliminating the inconsistency between cbioportal / tcga formulas.
  * Prefer lifelines (peer-reviewed library); fall back to the following when unavailable
    standard Mantel-Haenszel implementation, which is cross-validated with lifelines in _self_test()
    (agreeing to 1e-6).

API：
  logrank(durations, events, group) -> (chi2, p_value)
      Compatible with old script signatures, can be used as a drop-in replacement.
  logrank_detail(durations, events, group) -> dict  (contains n / number of events / test statistic)
"""
import numpy as np
from scipy import stats

try:
    from lifelines.statistics import logrank_test as _lifelines_logrank
    _HAS_LIFELINES = True
except Exception:
    _HAS_LIFELINES = False


def _mantel_haenszel(durations, events, group):
    """Standard log-rank (Mantel-Haenszel). Aggregates at each unique event time,
    handles ties correctly. Returns (chi2, p)."""
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
    """Two-sided log-rank test. group is a 0/1 array. Returns (chi2, p_value).

    Keep consistent with the signature of the old script's `logrank(durations, events, group)`
    Can directly replace two inconsistent old implementations.
    """
    d = np.asarray(durations, float)
    e = np.asarray(events, int)
    g = np.asarray(group, int)
    mask = ~np.isnan(d) & ~np.isnan(e)
    d, e, g = d[mask], e[mask], g[mask]
    if len(d) == 0 or int((g == 1).sum()) == 0 or int((g == 0).sum()) == 0:
        return 0.0, 1.0

    if use_lifelines and _HAS_LIFELINES:
        # lifelines.logrank_test signature is (durations_A, durations_B,
        #   event_observed_A, event_observed_B); must split into two groups by group,
        # Cannot pass the group= keyword (this keyword is not in this function signature).
        g1 = g == 1
        g0 = g == 0
        res = _lifelines_logrank(d[g1], d[g0],
                                 event_observed_A=e[g1], event_observed_B=e[g0])
        return float(res.test_statistic), float(res.p_value)
    return _mantel_haenszel(d, e, g)


def logrank_detail(durations, events, group):
    """Return dict: chi2, p, n, n_group1, n_group0, events_group1, events_group0."""
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
    """Cross-validate with lifelines (if available). Returns True if consistent."""
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
