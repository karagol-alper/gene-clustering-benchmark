"""
depbench.power
==============
What a clustered cohort can detect within groups.

Within-group permutation holds every group's class counts fixed.  For a score
with no ties the permutation distribution of the within-group concordance
AUC_W is then determined by those counts alone and not by the score, so the
test is distribution-free given the cohort structure and its null can be
computed once and reused.  That makes an exact power calculation cheap: draw
scores with a known population within-group AUC in the cohort's own structure
and count how often they fall outside the stored critical values.

Although the variant-effect application speaks of genes, nothing here is
specific to genes: `groups` may be any clustering of the rows.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

from depbench.metrics import auc_decompose

__all__ = ["within_null", "power_at", "power_curve", "min_detectable_within_auc",
           "delta_for_within_auc"]


def _within(y, s, g):
    return auc_decompose(y, s, g)["auc_within"]


def delta_for_within_auc(tau):
    """Mean shift that gives a population within-group AUC of tau when scores
    are standard normal: P(X_pos > X_neg) = Phi(delta / sqrt(2))."""
    return float(np.sqrt(2.0) * norm.ppf(tau))


def within_null(y, groups, n_perm=20000, seed=11):
    """Null distribution of AUC_W under within-group label permutation.

    Distribution-free for continuous scores, so any continuous score serves to
    generate it; `power_at` relies on that.  Returns the array of null values."""
    y = np.asarray(y)
    groups = np.asarray(groups)
    rng = np.random.default_rng(seed)
    s = rng.normal(size=len(y))
    mixed = [np.where(groups == g)[0] for g in np.unique(groups)
             if 0 < y[groups == g].sum() < (groups == g).sum()]
    out = np.empty(n_perm)
    for b in range(n_perm):
        yp = y.copy()
        for ix in mixed:
            yp[ix] = rng.permutation(y[ix])
        out[b] = _within(yp, s, groups)
    return out


def power_at(tau, y, groups, null=None, reps=4000, seed=11, alpha=0.05):
    """Power of the two-sided within-group permutation test against a score
    whose population within-group AUC is tau, in this cohort's exact group and
    class structure."""
    y = np.asarray(y)
    groups = np.asarray(groups)
    if null is None:
        null = within_null(y, groups, seed=seed)
    lo, hi = np.percentile(null, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    rng = np.random.default_rng(seed)
    delta = delta_for_within_auc(tau)
    hits = 0
    for _ in range(reps):
        s = rng.normal(size=len(y)) + delta * y
        w = _within(y, s, groups)
        hits += (w <= lo) or (w >= hi)
    return hits / reps


def power_curve(y, groups, taus, null=None, reps=4000, seed=11, alpha=0.05):
    """Power at each tau.  Returns a list of (tau, power) pairs."""
    if null is None:
        null = within_null(y, groups, seed=seed)
    return [(t, power_at(t, y, groups, null, reps, seed + i, alpha))
            for i, t in enumerate(taus)]


def min_detectable_within_auc(curve, target=0.8):
    """Smallest within-group AUC reaching `target` power, by linear
    interpolation on the power curve.  Returns None if the curve never does."""
    a = np.array([t for t, _ in curve], dtype=float)
    p = np.array([v for _, v in curve], dtype=float)
    above = np.where(p >= target)[0]
    if not len(above):
        return None
    i = int(above[0])
    if i == 0:
        return float(a[0])
    return float(np.interp(target, [p[i - 1], p[i]], [a[i - 1], a[i]]))
