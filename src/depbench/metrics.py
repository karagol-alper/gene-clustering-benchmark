"""
depbench.metrics
================
Discrimination statistics for gene-clustered variant cohorts.

* fast_auc            rank-based AUC (Mann-Whitney), tie-aware
* auc_decompose       exact split of the AUC into within-gene and between-gene
                      pair components (AUC = w*AUC_within + (1-w)*AUC_between)
* delong_test         paired DeLong test, variants treated as independent
* obuchowski_test     paired test with Obuchowski's clustered variance (t on I - 1 df), genes
                      as the independent units
* cluster_bootstrap   percentile intervals from resampling whole genes
* within_gene_perm_p  stratified (within-gene) permutation test of AUC_within
"""
from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.stats import rankdata


# ────────────────────────────────────────────────────────────────── AUC
def fast_auc(y, s):
    """Mann-Whitney AUC with mid-rank ties; equals sklearn's roc_auc_score."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(s, dtype=float)
    n1 = y.sum()
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(s)
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def _pair_matrix(y, s):
    y = np.asarray(y).astype(bool)
    s = np.asarray(s, dtype=float)
    sp, sn = s[y], s[~y]
    return (sp[:, None] > sn[None, :]) + 0.5 * (sp[:, None] == sn[None, :])


def auc_decompose(y, s, g):
    """
    Split the AUC into its within-gene and between-gene pair components.

    Every positive-negative pair is either within a gene (both variants from
    the same gene) or between genes.  The overall AUC is the pair-count
    weighted mean of the two components:
        AUC = w * AUC_within + (1 - w) * AUC_between,
    where w is the fraction of positive-negative pairs that are within-gene.
    AUC_within is the empirical covariate-adjusted AUC (Janes & Pepe 2008)
    with gene as the stratifying covariate; a score carrying only gene-level
    information has AUC_within = 0.5 exactly.
    """
    y = np.asarray(y).astype(bool)
    g = np.asarray(g)
    M = _pair_matrix(y, s)
    same = g[y][:, None] == g[~y][None, :]
    n_w, n_b = same.sum(), (~same).sum()
    auc_w = M[same].mean() if n_w else np.nan
    auc_b = M[~same].mean() if n_b else np.nan
    return dict(auc=M.mean(), auc_within=auc_w, auc_between=auc_b,
                n_pairs=int(M.size), n_within=int(n_w), n_between=int(n_b),
                w_within=n_w / M.size)


def per_gene_auc(y, s, g):
    """AUC within each gene that contains both classes."""
    y = np.asarray(y)
    s = np.asarray(s, dtype=float)
    g = np.asarray(g)
    res = {}
    for gi in np.unique(g):
        m = g == gi
        if 0 < y[m].sum() < m.sum():
            res[gi] = fast_auc(y[m], s[m])
    return res


# ─────────────────────────────────────────────── paired AUC comparisons
def _placements(y, s):
    """DeLong structural components: placement value of each positive among
    negatives (V10) and of each negative among positives (V01)."""
    M = _pair_matrix(y, s)
    return M.mean(axis=1), M.mean(axis=0)  # V10 (len m), V01 (len n)


def delong_test(y, p1, p2):
    """Two-sided DeLong test for two correlated AUCs (independent variants).
    Returns (auc1, auc2, z, p)."""
    y = np.asarray(y).astype(bool)
    m, n = y.sum(), (~y).sum()
    a1, b1 = _placements(y, p1)
    a2, b2 = _placements(y, p2)
    S = np.cov(np.vstack([a1, a2])) / m + np.cov(np.vstack([b1, b2])) / n
    var = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    d = a1.mean() - a2.mean()
    z = d / np.sqrt(var) if var > 0 else 0.0
    return a1.mean(), a2.mean(), z, 2 * stats.norm.sf(abs(z))


def obuchowski_var(y, scores, g):
    """
    Clustered covariance matrix of AUCs for one or more scores
    (Obuchowski 1997, Biometrics 53:567).  Genes are the independent units;
    variants within a gene may be arbitrarily correlated.
    Returns (aucs, covariance matrix).
    """
    y = np.asarray(y).astype(bool)
    g = np.asarray(g)
    genes = np.unique(g)
    I = len(genes)
    M, N = y.sum(), (~y).sum()
    gp, gn = g[y], g[~y]
    V10s, V01s, aucs = [], [], []
    for s in scores:
        v10, v01 = _placements(y, s)            # per-variant placements
        aucs.append(v10.mean())
        V10s.append(np.array([v10[gp == gi].sum() for gi in genes]))
        V01s.append(np.array([v01[gn == gi].sum() for gi in genes]))
    m_i = np.array([(gp == gi).sum() for gi in genes])
    n_i = np.array([(gn == gi).sum() for gi in genes])
    k = len(scores)
    C = np.zeros((k, k))
    for a in range(k):
        for b in range(k):
            ea = V10s[a] - m_i * aucs[a]
            eb = V10s[b] - m_i * aucs[b]
            fa = V01s[a] - n_i * aucs[a]
            fb = V01s[b] - n_i * aucs[b]
            S10 = I / ((I - 1) * M) * np.sum(ea * eb)
            S01 = I / ((I - 1) * N) * np.sum(fa * fb)
            S11 = I / (I - 1) * np.sum(ea * fb)
            S11r = I / (I - 1) * np.sum(eb * fa)
            C[a, b] = S10 / M + S01 / N + (S11 + S11r) / (M * N)
    return np.array(aucs), C


def obuchowski_test(y, p1, p2, g, small_sample=True):
    """Two-sided paired comparison of two AUCs with the clustered variance.
    With small_sample=True (default) the statistic is referred to a t
    distribution on I - 1 degrees of freedom, I being the number of genes,
    rather than to the standard normal; with few clusters the variance is
    itself estimated imprecisely and the normal reference is anticonservative.
    Returns (auc1, auc2, statistic, p, se_difference)."""
    aucs, C = obuchowski_var(y, [p1, p2], g)
    var = C[0, 0] + C[1, 1] - 2 * C[0, 1]
    d = aucs[0] - aucs[1]
    z = d / np.sqrt(var) if var > 0 else 0.0
    I = len(np.unique(np.asarray(g)))
    p = 2 * stats.t.sf(abs(z), I - 1) if small_sample else 2 * stats.norm.sf(abs(z))
    return aucs[0], aucs[1], z, p, np.sqrt(max(var, 0))


def delong_var_single(y, s):
    y = np.asarray(y).astype(bool)
    a, b = _placements(y, s)
    return a.mean(), a.var(ddof=1) / y.sum() + b.var(ddof=1) / (~y).sum()


# ─────────────────────────────────────────────────────────── bootstrap
def cluster_bootstrap(y, g, stat_fn, B=5000, seed=11, cluster=True):
    """
    Bootstrap distribution of stat_fn(y_b, idx_b, g_b).

    cluster=True resamples whole genes with replacement (Field & Welsh 2007);
    each draw receives its own cluster label, so two draws of the same gene
    are treated as distinct clusters.  cluster=False resamples variants
    (naive bootstrap).  Resamples containing one class are discarded.
    stat_fn receives (y_b, idx_b, g_b) and returns a float or 1-D array.
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    g = np.asarray(g)
    genes = np.unique(g)
    members = [np.where(g == gi)[0] for gi in genes]
    outs = []
    for _ in range(B):
        if cluster:
            pick = rng.integers(0, len(genes), len(genes))
            idx = np.concatenate([members[i] for i in pick])
            gb = np.concatenate([np.full(len(members[i]), j) for j, i in enumerate(pick)])
        else:
            idx = rng.integers(0, len(y), len(y))
            gb = g[idx]
        yb = y[idx]
        if yb.min() == yb.max():
            continue
        outs.append(stat_fn(yb, idx, gb))
    return np.array(outs, dtype=float)


def percentile_ci(a, level=0.95):
    a = np.asarray(a, dtype=float)
    a = a[~np.isnan(a)] if a.ndim == 1 else a
    lo, hi = (1 - level) / 2 * 100, (1 + level) / 2 * 100
    return np.nanpercentile(a, lo, axis=0), np.nanpercentile(a, hi, axis=0)


# ──────────────────────────────────────────── within-gene permutation
def within_gene_perm_p(y, s, g, n_perm=10000, seed=11):
    """
    Stratified (within-gene) permutation test of within-gene discrimination
    (Monte-Carlo: labels are shuffled within genes for n_perm replicates).
    Labels are shuffled within genes (preserving every gene's class counts),
    so the reference distribution holds gene-level prevalence fixed.
    Returns (AUC_within observed, one-sided p for AUC_within >= observed,
    two-sided p, null mean).
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(int)
    g = np.asarray(g)
    obs = auc_decompose(y, s, g)["auc_within"]
    mixed = [np.where(g == gi)[0] for gi in np.unique(g)
             if 0 < y[g == gi].sum() < (g == gi).sum()]
    null = np.empty(n_perm)
    yp = y.copy()
    for b in range(n_perm):
        for idx in mixed:
            yp[idx] = rng.permutation(y[idx])
        null[b] = auc_decompose(yp, s, g)["auc_within"]
    p1 = (np.sum(null >= obs - 1e-12) + 1) / (n_perm + 1)
    p2 = (np.sum(np.abs(null - 0.5) >= abs(obs - 0.5) - 1e-12) + 1) / (n_perm + 1)
    return obs, p1, p2, null.mean()


def perm_p(null, obs, greater=True):
    null = np.asarray(null)
    if greater:
        return (np.sum(null >= obs) + 1) / (len(null) + 1)
    return (np.sum(null <= obs) + 1) / (len(null) + 1)


def mc_se(p, n):
    """Monte-Carlo standard error of a permutation p-value."""
    return np.sqrt(p * (1 - p) / n)
