#!/usr/bin/env python3
"""
Step 7 - design requirements by simulation.

How many genes are needed before the AUC of a fixed, gene-agnostic score is
estimated to a given precision?  Cohorts of G genes with m = 5 variants each
are drawn from the random-intercept model of the simulation step (gene SD of
log-odds sigma_g in {1, 2, 4}; within-gene effect chosen so that the
population within-gene AUC is 0.65 or 0.75).  For each design the sampling
SD of two estimators is computed over replicates:
    pooled AUC          all positive-negative pairs
    within-gene AUC     pairs within the same gene only
and the 95% half-width (1.96 x SD) is fitted by least squares on a log-log
scale to find the number of genes required for +/-0.05 and +/-0.03.
"""
import itertools
import os

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from depbench.common import N_JOBS, SEED, out

GS = [20, 40, 80, 160, 320, 640, 1280]
SIGMAS = [1.0, 2.0, 4.0]
TARGETS = {0.65: None, 0.75: None}
M = 5
R = int(os.environ.get("DEPBENCH_DESIGNREPS", 300))
BASE = 0.41


def within_and_pooled(y, s, g):
    df = pd.DataFrame({"y": y, "s": s, "g": g})
    n1 = df.groupby("g").y.transform("sum")
    n = df.groupby("g").y.transform("size")
    mixed = (n1 > 0) & (n1 < n)
    d = df[mixed]
    if len(d) == 0:
        w = np.nan
    else:
        r = d.groupby("g").s.rank()
        U = (r[d.y == 1].groupby(d.g[d.y == 1]).sum() -
             d[d.y == 1].groupby("g").size() * (d[d.y == 1].groupby("g").size() + 1) / 2)
        pairs = d[d.y == 1].groupby("g").size() * d[d.y == 0].groupby("g").size()
        w = U.sum() / pairs.sum()
    rk = pd.Series(s).rank().values
    N1 = y.sum(); N0 = len(y) - N1
    pooled = (rk[y == 1].sum() - N1 * (N1 + 1) / 2) / (N1 * N0)
    return pooled, w, int(mixed.groupby(df.g).first().sum()) if len(d) else 0


def beta_for(target, rng):
    """Within-gene effect giving the target population within-gene AUC,
    found on a large within-gene sample (intercept 0, no gene effect)."""
    x = rng.normal(size=200000)
    lo, hi = 0.0, 5.0
    for _ in range(40):
        b = (lo + hi) / 2
        yy = rng.random(len(x)) < 1 / (1 + np.exp(-b * x))
        r = pd.Series(x).rank().values
        n1 = yy.sum(); n0 = len(yy) - n1
        auc = (r[yy].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
        lo, hi = (b, hi) if auc < target else (lo, b)
    return (lo + hi) / 2


def intercept(sigma, beta, rng):
    z = rng.normal(size=200000) * np.sqrt(sigma ** 2 + beta ** 2)
    lo, hi = -15.0, 15.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if (1 / (1 + np.exp(-(mid + z)))).mean() < BASE else (lo, mid)
    return (lo + hi) / 2


def run(sigma, target, G, seed):
    rng = np.random.default_rng(seed)
    b = beta_for(target, rng)
    a = intercept(sigma, b, rng)
    pooled, within, nmix = [], [], []
    for _ in range(R):
        u = rng.normal(0, sigma, G)
        g = np.repeat(np.arange(G), M)
        x = rng.normal(size=G * M)
        y = (rng.random(G * M) < 1 / (1 + np.exp(-(a + u[g] + b * x)))).astype(int)
        if y.min() == y.max():
            continue
        p, w, k = within_and_pooled(y, x, g)
        pooled.append(p); within.append(w); nmix.append(k)
    return dict(gene_sd=sigma, target_within_auc=target, n_genes=G, n_variants=G * M,
                pooled_mean=np.mean(pooled), pooled_sd=np.std(pooled, ddof=1),
                within_mean=np.nanmean(within), within_sd=np.nanstd(within, ddof=1),
                mixed_genes_mean=np.mean(nmix))


if __name__ == "__main__":
    cells = list(itertools.product(SIGMAS, TARGETS, GS))
    rows = Parallel(n_jobs=N_JOBS)(delayed(run)(s, t, G, SEED + i) for i, (s, t, G) in enumerate(cells))
    sim = pd.DataFrame(rows)
    sim.to_csv(out("design_simulation.csv"), index=False)
    req = []
    for (s, t), d in sim.groupby(["gene_sd", "target_within_auc"]):
        for est in ["pooled", "within"]:
            hw = 1.96 * d[f"{est}_sd"].values
            lg, lh = np.log(d.n_genes.values), np.log(hw)
            k, c = np.polyfit(lg, lh, 1)          # log hw = c + k log G
            for h in (0.05, 0.03):
                G_req = float(np.exp((np.log(h) - c) / k))
                req.append(dict(gene_sd=s, target_within_auc=t, estimator=est, half_width=h,
                                genes_required=int(np.ceil(G_req)),
                                variants_required=int(np.ceil(G_req) * M), slope=round(k, 3)))
    req = pd.DataFrame(req)
    req.to_csv(out("design_requirements.csv"), index=False)
    print(sim.round(3).to_string(index=False))
    print(req.to_string(index=False))
