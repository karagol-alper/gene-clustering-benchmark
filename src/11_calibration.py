#!/usr/bin/env python3
"""
Step 11 - calibration of fold-specific models, and estimators of between-gene
discrimination that do not compare predictions across folds.

Under leave-one-gene-out validation every between-gene pathogenic-benign pair
compares predictions made by two different fitted models.  Those models are
trained on different gene sets and therefore on different class balances, so
their predictions are on different scales.  This step asks how much of the
below-chance pooled gene-held-out AUC that fact explains.

  A  calibration of the out-of-fold predictions under each design: intercept
     and slope of the logistic recalibration of the labels on the predicted
     log-odds, Brier score, and the relation between each fold's mean
     prediction and its training prevalence
  B  two fold-wise corrections, each a monotone transformation applied
     separately within each fold, so that neither can change the within-gene
     component:
       prior shift     logit(p) - logit(pi_train) + logit(pi_cohort)
       fold rank       percentile of each test prediction among the training
                       predictions of the same fold
  C  leave-two-genes-out: for every unordered pair of genes, a model is fitted
     on the remaining genes and used to score both held-out genes, so that
     each between-gene pair is scored by a single model that saw neither of
     its genes.  This is the estimator of between-gene discrimination for new
     genes that is free of cross-fold comparison.
"""
import json

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from depbench.calibration import (calibration_fit, fold_rank,
                                  leave_two_groups_out_auc, prior_shift)
from depbench.common import (ALL_COLS, DEP_COLS, N_BOOT, N_JOBS, SEED,
                             cross_validate, gene_prevalence,
                             gene_prevalence_logo, make_folds, out)
from depbench.metrics import (auc_decompose, cluster_bootstrap, fast_auc,
                              percentile_ci)

lab = pd.read_csv(out("cohort_full.csv"))
y, g = lab.label.values, lab.Gene_Name.values
genes = np.unique(g)
PI = y.mean()
EPS = 1e-6

MODELS = {
    "Position-level scores": (DEP_COLS, None),
    "Full 19-feature model": (ALL_COLS, 9),
    "AlphaMissense only": (["AM_Specific_Score"], None),
}


# ── A. calibration of out-of-fold predictions ───────────────────────────
rows, fold_rows = [], []
preds_store = {}
for name, (cols, k) in MODELS.items():
    X = lab[cols].values.astype(float)
    for split in ("variant", "gene"):
        p, _, _ = cross_validate(X, y, g, split, "lr", k=k)
        d = auc_decompose(y, p, g)
        rows.append(dict(model=name, split=split, auc=round(d["auc"], 4),
                         auc_within=round(d["auc_within"], 4),
                         auc_between=round(d["auc_between"], 4),
                         **{kk: round(vv, 4) for kk, vv in calibration_fit(y, p).items()}))
        preds_store[(name, split)] = p
        print(rows[-1], flush=True)
        if split == "gene":
            for train, test in make_folds(y, g, "gene"):
                fold_rows.append(dict(model=name, gene=g[test][0], n=len(test),
                                      train_prevalence=round(float(y[train].mean()), 4),
                                      gene_prevalence=round(float(y[test].mean()), 4),
                                      mean_prediction=round(float(p[test].mean()), 4)))
cal = pd.DataFrame(rows)
cal.to_csv(out("calibration_by_design.csv"), index=False)
folds_df = pd.DataFrame(fold_rows)
folds_df.to_csv(out("fold_calibration.csv"), index=False)
res = {}
f0 = folds_df[folds_df.model == "Position-level scores"]
res["fold_mean_pred_vs_train_prevalence_r"] = float(np.corrcoef(f0.train_prevalence,
                                                                f0.mean_prediction)[0, 1])
res["fold_mean_pred_vs_gene_prevalence_r"] = float(np.corrcoef(f0.gene_prevalence,
                                                               f0.mean_prediction)[0, 1])
res["train_vs_gene_prevalence_r"] = float(np.corrcoef(f0.gene_prevalence,
                                                      f0.train_prevalence)[0, 1])
print("fold correlations:", res, flush=True)


# ── B. fold-wise monotone corrections ───────────────────────────────────
def training_predictions(X, k):
    """In-fold predictions on each fold's own training set."""
    outp = []
    for train, test in make_folds(y, g, "gene"):
        pt, _, _ = cross_validate(X[train], y[train], g[train], "gene", "lr", k=k)
        outp.append(pt)
    return outp


corr_rows = []
folds_g = make_folds(y, g, "gene")
for name, (cols, k) in MODELS.items():
    X = lab[cols].values.astype(float)
    p = preds_store[(name, "gene")]
    tp = training_predictions(X, k)
    for label, q in [("none", p), ("prior shift", prior_shift(p, folds_g, y)),
                     ("fold rank", fold_rank(p, folds_g, tp))]:
        d = auc_decompose(y, q, g)
        lo, hi = percentile_ci(cluster_bootstrap(
            y, g, lambda yb, idx, gb, q=q: fast_auc(yb, q[idx]), B=N_BOOT, seed=SEED))
        corr_rows.append(dict(model=name, correction=label, auc=round(d["auc"], 4),
                              auc_within=round(d["auc_within"], 4),
                              auc_between=round(d["auc_between"], 4),
                              cluster_ci=f"{lo:.3f}-{hi:.3f}"))
        print(corr_rows[-1], flush=True)

# positive control: a predictor that knows only the training prevalence
prev_logo = gene_prevalence_logo(lab)
q = prior_shift(np.clip(prev_logo, EPS, 1 - EPS), folds_g, y)
corr_rows.append(dict(model="Training prevalence only", correction="none",
                      auc=round(fast_auc(y, prev_logo), 4),
                      auc_within=round(auc_decompose(y, prev_logo, g)["auc_within"], 4),
                      auc_between=round(auc_decompose(y, prev_logo, g)["auc_between"], 4),
                      cluster_ci="—"))
corr_rows.append(dict(model="Training prevalence only", correction="prior shift",
                      auc=round(fast_auc(y, q), 4),
                      auc_within=round(auc_decompose(y, q, g)["auc_within"], 4),
                      auc_between=round(auc_decompose(y, q, g)["auc_between"], 4),
                      cluster_ci="—"))
pd.DataFrame(corr_rows).to_csv(out("fold_corrections.csv"), index=False)


# ── C. leave-two-genes-out between-gene discrimination ──────────────────
def lpgo_pairs(name, cols, k):
    """Between-gene AUC with every pair ranked by one model blind to both of
    its genes (depbench.calibration.leave_two_groups_out_auc)."""
    X = lab[cols].values.astype(float)

    def fit_predict(train, test):
        return cross_validate(X, y, g, "gene", "lr", k=k, folds=[(train, test)])[0][test]

    auc, n = leave_two_groups_out_auc(fit_predict, X, y, g)
    return dict(model=name, auc_between_lpgo=round(auc, 4), n_pairs=n)


lp = Parallel(n_jobs=N_JOBS)(delayed(lpgo_pairs)(n, c, k) for n, (c, k) in MODELS.items())
lp = pd.DataFrame(lp)
logo_b = {n: auc_decompose(y, preds_store[(n, "gene")], g)["auc_between"] for n in MODELS}
lp["auc_between_logo"] = [round(logo_b[n], 4) for n in lp.model]
lp["difference"] = (lp.auc_between_lpgo - lp.auc_between_logo).round(4)
lp.to_csv(out("leave_two_genes_out.csv"), index=False)
print(lp.to_string(index=False))

res["n_gene_pairs"] = int(len(list(itertools.combinations(genes, 2))))
json.dump(res, open(out("calibration.json"), "w"), indent=1, default=float)
print("calibration step written")
