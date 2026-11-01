#!/usr/bin/env python3
"""
Step 12 - how far the conclusions depend on single genes, what this cohort
can detect, and direct comparisons between scores within genes.

  A  gene-drop influence: every analysis that produces a pooled gene-held-out
     AUC or a within-gene AUC is repeated with each gene removed from the
     cohort, so that the dependence on any one gene is explicit
  B  detectable within-gene effect: because the within-gene permutation null
     of AUC_W depends only on each gene's class counts and not on the score
     (the test is distribution-free for continuous scores), the null is
     computed once and used to obtain the power of the stratified permutation
     test against scores of known within-gene AUC in this cohort's exact
     structure, and the smallest within-gene AUC detectable at 80% power
  C  paired comparisons of within-gene discrimination between scores, with
     stratified bootstrap intervals for the difference
"""
import json

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from depbench.common import DEP_COLS, N_BOOT, N_JOBS, SEED, cross_validate, out
from depbench.metrics import auc_decompose, fast_auc, percentile_ci
from depbench.power import (delta_for_within_auc, min_detectable_within_auc,
                            power_at, within_null)

lab = pd.read_csv(out("cohort_full.csv"))
y, g = lab.label.values, lab.Gene_Name.values
genes = np.unique(g)
mixed = [x for x in genes if 0 < y[g == x].sum() < (g == x).sum()]
rng = np.random.default_rng(SEED)
res = {}

SCORES = {"AlphaMissense (variant)": "AM_Specific_Score", "DCS": "DCS_Score",
          "RCS": "RCS_Score", "Conservation": "conservation",
          "Normalised mobility": "Norm_Dynamics"}


def within_auc(yv, s, gv):
    return auc_decompose(yv, s, gv)["auc_within"]


# ── A. gene-drop influence ──────────────────────────────────────────────
def drop_gene(gn):
    keep = g != gn
    yk, gk = y[keep], g[keep]
    row = dict(dropped=gn, n_removed=int((~keep).sum()),
               n_pathogenic_removed=int(y[~keep].sum()))
    if len(np.unique(yk)) == 2 and len(np.unique(gk)) > 1:
        p, _, _ = cross_validate(lab.loc[keep, DEP_COLS].values, yk, gk, "gene", "lr", k=None)
        d = auc_decompose(yk, p, gk)
        row.update(pos_auc=round(d["auc"], 4), pos_within=round(d["auc_within"], 4),
                   pos_between=round(d["auc_between"], 4))
    for name, col in SCORES.items():
        s = lab.loc[keep, col].values.astype(float)
        w = within_auc(yk, s, gk)
        row[f"within_{name}"] = round(w, 4) if np.isfinite(w) else np.nan
    d0 = auc_decompose(yk, lab.loc[keep, "AM_Specific_Score"].values, gk)
    row["n_within_pairs"] = int(d0["n_within"])
    return row


infl = pd.DataFrame(Parallel(n_jobs=N_JOBS)(delayed(drop_gene)(gn) for gn in genes))
infl = infl.sort_values("n_within_pairs")
infl.to_csv(out("gene_drop_influence.csv"), index=False)
print(infl.to_string(index=False), flush=True)
full = {name: within_auc(y, lab[col].values.astype(float), g) for name, col in SCORES.items()}
res["within_full_cohort"] = {k: float(v) for k, v in full.items()}
res["within_drop_PGK1"] = {k: (None if pd.isna(v) else float(v)) for k, v in
                           infl[infl.dropped == "PGK1"].iloc[0][
                               [f"within_{n}" for n in SCORES]].items()}
res["n_within_pairs_drop_PGK1"] = int(infl[infl.dropped == "PGK1"].n_within_pairs.iloc[0])
res["pos_auc_range_over_drops"] = [float(infl.pos_auc.min()), float(infl.pos_auc.max())]
res["pos_auc_most_influential"] = str(infl.loc[infl.pos_auc.idxmax(), "dropped"])


# ── B. detectable within-gene effect ────────────────────────────────────
null = within_null(y, g, 20000, SEED)
# the null is distribution-free: an independent run must give the same quantiles
null2 = within_null(y, g, 20000, SEED + 1)
q_lo, q_hi = np.percentile(null, [2.5, 97.5])
res["null_check_max_quantile_difference"] = float(
    np.max(np.abs(np.percentile(null, [1, 2.5, 5, 50, 95, 97.5, 99]) -
                  np.percentile(null2, [1, 2.5, 5, 50, 95, 97.5, 99]))))
res["within_null_mean"] = float(null.mean())
res["within_null_95"] = [float(q_lo), float(q_hi)]


TAUS = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]
pw = Parallel(n_jobs=N_JOBS)(delayed(power_at)(t, y, g, null, 4000, SEED + i)
                            for i, t in enumerate(TAUS))
power = pd.DataFrame(dict(target_within_auc=TAUS, power=[round(x, 4) for x in pw]))
power.to_csv(out("within_gene_power.csv"), index=False)
print(power.to_string(index=False), flush=True)
a = np.array(TAUS)
pv = np.array(pw)
res["min_detectable_within_auc_80"] = min_detectable_within_auc(list(zip(TAUS, pw)), 0.8)
res["power_at_065"] = float(pv[a == 0.65][0])
res["power_at_070"] = float(pv[a == 0.70][0])


# ── C. paired within-gene comparisons ───────────────────────────────────
strata = [np.where((g == x) & (y == c))[0] for x in mixed for c in (0, 1)]


def strat_boot_diff(s1, s2, B=N_BOOT, seed=SEED):
    r = np.random.default_rng(seed)
    d = np.empty(B)
    for b in range(B):
        idx = np.concatenate([r.choice(ix, len(ix)) for ix in strata])
        d[b] = within_auc(y[idx], s1[idx], g[idx]) - within_auc(y[idx], s2[idx], g[idx])
    return d


PAIRS = [("DCS", "Conservation"), ("RCS", "Conservation"), ("DCS", "RCS"),
         ("Normalised mobility", "Conservation"),
         ("AlphaMissense (variant)", "DCS"), ("AlphaMissense (variant)", "Conservation")]
rows = []
for a_, b_ in PAIRS:
    s1 = lab[SCORES[a_]].values.astype(float)
    s2 = lab[SCORES[b_]].values.astype(float)
    w1, w2 = within_auc(y, s1, g), within_auc(y, s2, g)
    bs = strat_boot_diff(s1, s2)
    lo, hi = percentile_ci(bs)
    rows.append(dict(comparison=f"{a_} vs {b_}", within_1=round(w1, 4), within_2=round(w2, 4),
                     difference=round(w1 - w2, 4), ci_lo=round(float(lo), 4),
                     ci_hi=round(float(hi), 4), ci_width=round(float(hi - lo), 4),
                     excludes_zero=bool(lo > 0 or hi < 0)))
    print(rows[-1], flush=True)
pw2 = pd.DataFrame(rows)
pw2.to_csv(out("within_gene_paired.csv"), index=False)
res["median_paired_ci_width"] = float(pw2.ci_width.median())
res["n_paired_excluding_zero"] = int(pw2.excludes_zero.sum())

json.dump(res, open(out("robustness.json"), "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
