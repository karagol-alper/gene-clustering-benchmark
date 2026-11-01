"""
depbench.fastlr
===============
Vectorised L2-penalised logistic regression with balanced class weights,
solving the same objective as
    sklearn.linear_model.LogisticRegression(C=C, class_weight="balanced")
(0.5*||w||^2 + C * sum_i c_i * logloss_i, intercept unpenalised) by Newton's
method.  Used where many thousands of small fits are required (simulations
and permutation nulls); agreement with scikit-learn is checked in tests.
"""
import numpy as np


def fit_lr(X, y, C=1.0, balanced=True, n_iter=50, tol=1e-10):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n, p = X.shape
    if balanced:
        n1 = y.sum(); n0 = n - n1
        cw = np.where(y == 1, n / (2 * n1), n / (2 * n0))
    else:
        cw = np.ones(n)
    Z = np.hstack([np.ones((n, 1)), X])
    beta = np.zeros(p + 1)
    P = np.eye(p + 1) / C
    P[0, 0] = 0.0
    for _ in range(n_iter):
        eta = Z @ beta
        mu = 1 / (1 + np.exp(-eta))
        g = Z.T @ (cw * (mu - y)) + P @ beta
        W = cw * mu * (1 - mu)
        H = (Z * W[:, None]).T @ Z + P
        step = np.linalg.solve(H + 1e-12 * np.eye(p + 1), g)
        beta -= step
        if np.max(np.abs(step)) < tol:
            break
    return beta


def predict_lr(beta, X):
    X = np.atleast_2d(np.asarray(X, dtype=float))
    return 1 / (1 + np.exp(-(beta[0] + X @ beta[1:])))


def standardise(Xtr, Xte):
    mu = Xtr.mean(0); sd = Xtr.std(0)
    sd[sd == 0] = 1.0
    return (Xtr - mu) / sd, (Xte - mu) / sd


def cv_predict(X, y, folds, C=1.0, scale=True):
    """Out-of-fold predictions for arbitrary folds (list of (train, test))."""
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    pred = np.full(len(y), np.nan)
    for tr, te in folds:
        if y[tr].min() == y[tr].max():
            pred[te] = y[tr].mean()
            continue
        Xtr, Xte = (standardise(X[tr], X[te]) if scale else (X[tr], X[te]))
        b = fit_lr(Xtr, y[tr], C=C)
        pred[te] = predict_lr(b, Xte)
    return pred


def loo_folds(n):
    idx = np.arange(n)
    return [(np.delete(idx, i), np.array([i])) for i in range(n)]


def logo_folds(groups):
    groups = np.asarray(groups)
    return [(np.where(groups != g)[0], np.where(groups == g)[0]) for g in np.unique(groups)]
