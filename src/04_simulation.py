#!/usr/bin/env python3
"""
Step 4 - simulation with known estimands.

Generative model (gene j, variant i):
    u_j  ~ N(0, sigma_g^2)          gene-level log-odds offset (prevalence)
    f_j  ~ N(0, 1)                  gene-level feature offset (fingerprint)
    x_ij ~ N(0, 1)                  variant-level property
    y_ij ~ Bernoulli(expit(a + u_j + beta * x_ij))
    observed feature  x~_ij = x_ij + gamma * f_j
The intercept a is set so that the marginal pathogenic rate is ~0.41.

Predictors evaluated in every replicate
    gene prevalence     the pathogenic fraction of the variant's gene
    feature model       logistic regression on x~ (balanced class weights)
    local model         7-nearest-neighbour vote on (x~, h~), where
                        h~_ij = f2_j + 0.3 e_ij is a second, purely gene-level
                        position descriptor; this learner can identify genes
                        from the two descriptors, as a random forest can

Estimators
    variant-LOO   pooled AUC of leave-one-variant-out predictions
    gene-LOO      pooled AUC of leave-one-gene-out predictions
    gene-LOO (w)  within-gene component of the leave-one-gene-out AUC

Estimands (model refitted on the complete simulated cohort and applied to
large independent test sets)
    seen-gene AUC    new variants drawn from the same genes
    new-gene AUC     variants from new genes
    within-gene AUC  within-gene pairs among new-gene variants

Outputs: simulation.csv (cell means/SDs), simulation_replicates.csv.gz
"""
import itertools
import os

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from depbench.common import N_JOBS, SEED, out
from depbench.fastlr import fit_lr, predict_lr
from depbench.metrics import auc_decompose, fast_auc

R = int(os.environ.get("DEPBENCH_SIMREPS", 150))
COHORT_SIZES = [34, 13, 9, 5, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1]
SIGMAS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
DESIGNS = {           # name: list of gene sizes
    "12 genes x 2": [2] * 12, "12 genes x 5": [5] * 12, "12 genes x 10": [10] * 12,
    "48 genes x 5": [5] * 48, "cohort-matched (18 genes)": COHORT_SIZES,
}
BETAS = [0.0, 0.75]
GAMMAS = [0.0, 1.0]
BASE = 0.41


def intercept_for(sigma, beta, rng):
    """Intercept giving a marginal pathogenic rate of BASE (Monte Carlo)."""
    z = rng.normal(size=200000) * np.sqrt(sigma ** 2 + beta ** 2)
    lo, hi = -15.0, 15.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if (1 / (1 + np.exp(-(mid + z)))).mean() < BASE:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def draw(sizes, sigma, beta, gamma, a, rng, u=None, f=None):
    G = len(sizes)
    u = rng.normal(0, sigma, G) if u is None else u
    f = rng.normal(0, 1, (G, 2)) if f is None else f
    g = np.repeat(np.arange(G), sizes)
    x = rng.normal(size=len(g))
    y = (rng.random(len(g)) < 1 / (1 + np.exp(-(a + u[g] + beta * x)))).astype(int)
    h = f[g, 1] + 0.3 * rng.normal(size=len(g))
    return g, x + gamma * f[g, 0], y, u, f, h


def knn_predict(Ztr, ytr, Zte, k=7, exclude=None):
    """Distance-weighted k-NN probability.  `exclude` is a boolean matrix
    (n_te x n_tr) of training points that may not be used for a test point."""
    mu, sd = Ztr.mean(0), Ztr.std(0)
    sd[sd == 0] = 1.0
    A, B = (Ztr - mu) / sd, (Zte - mu) / sd
    D = np.sqrt(((B[:, None, :] - A[None, :, :]) ** 2).sum(-1))
    if exclude is not None:
        D = np.where(exclude, np.inf, D)
    idx = np.argpartition(D, k, axis=1)[:, :k]
    d = np.take_along_axis(D, idx, 1)
    w = 1 / (d + 1e-6)
    return (w * ytr[idx]).sum(1) / w.sum(1)


def lr_fit_predict(xtr, ytr, xte):
    if ytr.min() == ytr.max():
        return np.full(len(xte), ytr.mean())
    mu, sd = xtr.mean(), xtr.std() or 1.0
    b = fit_lr(((xtr - mu) / sd)[:, None], ytr)
    return predict_lr(b, ((xte - mu) / sd)[:, None])


def one_rep(sizes, sigma, beta, gamma, a, seed):
    rng = np.random.default_rng(seed)
    for _ in range(100):                       # need both classes in the cohort
        g, xt, y, u, f, h = draw(sizes, sigma, beta, gamma, a, rng)
        if 0 < y.sum() < len(y):
            break
    n = len(y)
    res = {}
    # ---- feature model: CV estimators
    pv = np.empty(n)
    for i in range(n):
        tr = np.arange(n) != i
        pv[i] = lr_fit_predict(xt[tr], y[tr], xt[i:i + 1])[0]
    pg = np.empty(n)
    for gj in np.unique(g):
        te = g == gj
        pg[te] = lr_fit_predict(xt[~te], y[~te], xt[te])
    dg = auc_decompose(y, pg, g)
    res.update(feat_vloo=fast_auc(y, pv), feat_gloo=dg["auc"], feat_gloo_w=dg["auc_within"],
               feat_gloo_b=dg["auc_between"])
    # ---- local (k-NN) model on two descriptors: CV estimators
    Z = np.c_[xt, h]
    eye = np.eye(n, dtype=bool)
    kv = knn_predict(Z, y, Z, exclude=eye)                 # leave-one-variant-out
    kg = knn_predict(Z, y, Z, exclude=(g[:, None] == g[None, :]))   # leave-one-gene-out
    dk = auc_decompose(y, kg, g)
    res.update(knn_vloo=fast_auc(y, kv), knn_gloo=dk["auc"], knn_gloo_w=dk["auc_within"])
    # ---- gene-prevalence lookup: CV estimators
    prev_v = np.empty(n)
    for i in range(n):
        same = (g == g[i]); same[i] = False
        prev_v[i] = y[same].mean() if same.any() else y.mean()
    prev_g = np.array([y[g != gj].mean() for gj in g])
    res.update(prev_vloo=fast_auc(y, prev_v), prev_gloo=fast_auc(y, prev_g))
    # ---- estimands: refit on the complete cohort, apply to new data
    rng2 = np.random.default_rng(seed + 10_000_019)
    # seen genes: 40 new variants per existing gene
    gs, xs, ys, _, _, hs = draw([40] * len(sizes), sigma, beta, gamma, a, rng2, u=u, f=f)
    # new genes: 400 genes x 10 variants
    gn, xn, yn, _, _, hn = draw([10] * 400, sigma, beta, gamma, a, rng2)
    p_seen = lr_fit_predict(xt, y, xs)
    p_new = lr_fit_predict(xt, y, xn)
    dn = auc_decompose(yn, p_new, gn)
    res.update(feat_true_seen=fast_auc(ys, p_seen), feat_true_new=dn["auc"],
               feat_true_within=dn["auc_within"])
    ks = knn_predict(Z, y, np.c_[xs, hs])
    kn = knn_predict(Z, y, np.c_[xn, hn])
    dkn = auc_decompose(yn, kn, gn)
    res.update(knn_true_seen=fast_auc(ys, ks), knn_true_new=dkn["auc"],
               knn_true_within=dkn["auc_within"])
    prev_by_gene = np.array([y[g == j].mean() for j in range(len(sizes))])
    res.update(prev_true_seen=fast_auc(ys, prev_by_gene[gs]), prev_true_new=0.5)
    res["n_mixed_genes"] = int(sum(0 < y[g == j].sum() < (g == j).sum() for j in range(len(sizes))))
    return res


def run_cell(design, sigma, beta, gamma, cell_seed):
    rng = np.random.default_rng(cell_seed)
    a = intercept_for(sigma, beta, rng)
    sizes = DESIGNS[design]
    rows = []
    for r in range(R):
        d = one_rep(sizes, sigma, beta, gamma, a, cell_seed * 1000 + r)
        d.update(design=design, gene_sd=sigma, beta=beta, gamma=gamma, rep=r)
        rows.append(d)
    return rows


if __name__ == "__main__":
    cells = list(itertools.product(DESIGNS, SIGMAS, BETAS, GAMMAS))
    res = Parallel(n_jobs=N_JOBS, verbose=0)(
        delayed(run_cell)(d, s, b, gm, SEED * 7919 + i) for i, (d, s, b, gm) in enumerate(cells))
    reps = pd.DataFrame([r for cell in res for r in cell])
    reps.to_csv(out("simulation_replicates.csv.gz"), index=False)
    reps["feat_infl_vloo_new"] = reps.feat_vloo - reps.feat_true_new
    reps["feat_bias_vloo_seen"] = reps.feat_vloo - reps.feat_true_seen
    reps["feat_bias_gloo_new"] = reps.feat_gloo - reps.feat_true_new
    reps["feat_bias_gloo_w"] = reps.feat_gloo_w - reps.feat_true_within
    reps["feat_vloo_minus_gloo"] = reps.feat_vloo - reps.feat_gloo
    reps["knn_infl_vloo_new"] = reps.knn_vloo - reps.knn_true_new
    reps["knn_bias_vloo_seen"] = reps.knn_vloo - reps.knn_true_seen
    reps["knn_bias_gloo_new"] = reps.knn_gloo - reps.knn_true_new
    reps["knn_bias_gloo_w"] = reps.knn_gloo_w - reps.knn_true_within
    reps["prev_infl_vloo_new"] = reps.prev_vloo - reps.prev_true_new
    reps["prev_bias_gloo_new"] = reps.prev_gloo - reps.prev_true_new
    keys = ["design", "gene_sd", "beta", "gamma"]
    agg = reps.groupby(keys).agg(["mean", "std"])
    agg.columns = [f"{a}_{b}" for a, b in agg.columns]
    agg = agg.reset_index()
    agg.to_csv(out("simulation.csv"), index=False)
    show = ["design", "gene_sd", "beta", "gamma", "feat_vloo_mean", "feat_gloo_mean",
            "feat_true_seen_mean", "feat_true_new_mean", "feat_infl_vloo_new_mean",
            "feat_bias_gloo_new_mean", "feat_bias_gloo_w_mean", "knn_vloo_mean", "knn_gloo_mean",
            "knn_true_seen_mean", "knn_true_new_mean", "prev_vloo_mean", "prev_gloo_mean"]
    print(agg[show].round(3).to_string(index=False))
