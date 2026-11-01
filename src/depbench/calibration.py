"""
depbench.calibration
====================
Calibration of fold-specific models, and estimators of between-group
discrimination that do not compare predictions made by different models.

Under leave-one-group-out validation every between-group pair of a positive
and a negative compares predictions from two different fitted models.  Those
models are trained on different subsets and therefore on different class
balances, so their predictions sit on different scales and the pooled AUC
absorbs the difference.  The corrections here act separately within each fold
and are strictly increasing, so none of them can change the within-group
component of the AUC; `leave_two_groups_out_auc` removes the cross-fold
comparison altogether.

Although the variant-effect application speaks of genes, nothing here is
specific to genes: `groups` may be any clustering of the rows.
"""
from __future__ import annotations

import itertools

import numpy as np
from scipy.special import expit, logit
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression

EPS = 1e-6

__all__ = ["calibration_fit", "prior_shift", "fold_rank", "training_predictions",
           "leave_two_groups_out_auc"]


def _clip(p):
    return np.clip(np.asarray(p, dtype=float), EPS, 1 - EPS)


def calibration_fit(y, p):
    """Intercept and slope of logit P(y = 1) = a + b * logit(p), with the Brier
    score.  Perfect calibration is a = 0, b = 1; b < 1 marks an overconfident
    predictor and b < 0 a reversed one."""
    y = np.asarray(y, dtype=float)
    z = logit(_clip(p))[:, None]
    m = LogisticRegression(C=1e6, solver="lbfgs", max_iter=1000).fit(z, y)
    return dict(cal_intercept=float(m.intercept_[0]), cal_slope=float(m.coef_[0, 0]),
                brier=float(np.mean((_clip(p) - y) ** 2)), mean_pred=float(np.mean(p)))


def prior_shift(p, folds, y, pi=None):
    """Remove each fold's training prior from its predicted log-odds:

        logit(p') = logit(p) - logit(pi_train) + logit(pi)

    `folds` is a list of (train_index, test_index) pairs.  The map is strictly
    increasing within a fold, so the within-group component of the AUC is
    unchanged.  Applied to a predictor that scores every row by its fold's
    training prevalence, it returns an AUC of exactly 0.5."""
    y = np.asarray(y, dtype=float)
    p = _clip(p)
    pi = float(y.mean()) if pi is None else float(pi)
    q = p.copy()
    for train, test in folds:
        pit = float(np.clip(y[train].mean(), EPS, 1 - EPS))
        q[test] = expit(logit(p[test]) - logit(pit) + logit(pi))
    return q


def training_predictions(fit_predict, X, y, groups, folds):
    """In-fold predictions on each fold's own training rows, for `fold_rank`.

    `fit_predict(X_sub, y_sub, groups_sub)` must return out-of-fold predictions
    for the rows it is given."""
    return [fit_predict(X[train], y[train], groups[train]) for train, _ in folds]


def fold_rank(p, folds, train_preds):
    """Replace each test prediction by its percentile among the predictions the
    same model makes on its own training rows.  Ties created by the finite
    training set are broken by the prediction itself, so the map is strictly
    increasing within a fold and the within-group component is unchanged."""
    p = np.asarray(p, dtype=float)
    q = p.copy()
    for i, (_, test) in enumerate(folds):
        ref = np.sort(np.asarray(train_preds[i], dtype=float))
        pc = np.searchsorted(ref, p[test], side="left") / max(len(ref), 1)
        q[test] = pc + 1e-9 * (rankdata(p[test]) / (len(test) + 1.0))
    return q


def leave_two_groups_out_auc(fit_predict, X, y, groups):
    """Between-group AUC in which each pair is ranked by one model that saw
    neither of its groups.

    For every unordered pair of groups a model is fitted on the complement and
    used to score both held-out groups; concordance is then counted over
    positive-negative pairs whose members come from the two different held-out
    groups.  Summing over all pairs covers every between-group pair exactly
    once.  Costs one fit per pair of groups.

    `fit_predict(train_idx, test_idx)` must return predictions for `test_idx`.
    Returns (auc, n_pairs)."""
    y = np.asarray(y)
    groups = np.asarray(groups)
    uniq = np.unique(groups)
    num = den = 0.0
    for ga, gb in itertools.combinations(uniq, 2):
        held = (groups == ga) | (groups == gb)
        train, test = np.where(~held)[0], np.where(held)[0]
        if len(np.unique(y[train])) < 2:
            continue
        p = np.asarray(fit_predict(train, test), dtype=float)
        pos = test[y[test] == 1]
        neg = test[y[test] == 0]
        pm = {int(i): float(v) for i, v in zip(test, p)}
        for i in pos:
            for j in neg:
                if groups[i] == groups[j]:
                    continue
                a, b = pm[int(i)], pm[int(j)]
                num += 1.0 if a > b else (0.5 if a == b else 0.0)
                den += 1
    return (num / den if den else np.nan), int(den)
