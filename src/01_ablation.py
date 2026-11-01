#!/usr/bin/env python3
"""
Step 1 - factorial validation ablation.

Four binary factors, each a pipeline stage that can either see or not see the
held-out data:
    split      variant-level (leave-one-variant-out) or gene-level
               (leave-one-gene-out) resampling
    scaling    median imputation + Yeo-Johnson fitted on the complete labelled
               set ('full set') or inside each training fold ('in fold')
    selection  univariate ANOVA-F selection of nine features, full set or in fold
    tuning     regularisation strength C chosen by repeated stratified 5-fold
               AUC on the complete labelled set ('full set') or by an inner
               cross-validation nested inside each training fold ('in fold';
               inner folds grouped by gene whenever the outer split is by gene)

Logistic regression: full 2^4 design (16 configurations).
Random forest and gradient boosting: 2^3 design (split x scaling x selection)
with fixed hyperparameters, as in the main analysis.
Outputs: ablation_lr16.csv, ablation_full.csv, ablation_effects.csv,
         ablation_predictions.csv
"""
import itertools

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import (GroupKFold, LeaveOneGroupOut, LeaveOneOut,
                                     RepeatedStratifiedKFold, StratifiedKFold)
from sklearn.pipeline import Pipeline

from depbench.common import (ALL_COLS, N_JOBS, SEED, FoldCache, cross_validate,
                             out, preprocessor)
from depbench.metrics import auc_decompose

C_GRID = np.logspace(-3, 2, 11)
lab = pd.read_csv(out("cohort_full.csv"))
X, y, g = lab[ALL_COLS].values.astype(float), lab.label.values, lab.Gene_Name.values


def lr():
    return LogisticRegression(class_weight="balanced", max_iter=2000)


def inner_auc_all(tr, Xsrc, cache, leak_scale, leak_select, grouped):
    """Inner-CV AUC for every C in the grid (preprocessing cached per inner
    split; selection and model refitted)."""
    ytr, gtr = y[tr], g[tr]
    if grouped:
        inner = list(GroupKFold(n_splits=5).split(tr, ytr, gtr))
    else:
        inner = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(tr, ytr))
    scores = np.zeros(len(C_GRID))
    count = 0
    for a, b in inner:
        if ytr[b].min() == ytr[b].max():
            continue
        if leak_scale:
            Xa, Xb = Xsrc[tr[a]], Xsrc[tr[b]]
        else:
            Xa, Xb = cache.get(tr[a], tr[b])
        if not leak_select and Xa.shape[1] > 9:
            sel = SelectKBest(f_classif, k=9).fit(Xa, ytr[a])
            Xa, Xb = sel.transform(Xa), sel.transform(Xb)
        for j, C in enumerate(C_GRID):
            m = LogisticRegression(class_weight="balanced", max_iter=2000, C=C).fit(Xa, ytr[a])
            scores[j] += roc_auc_score(ytr[b], m.predict_proba(Xb)[:, 1])
        count += 1
    return scores / max(count, 1)


def leaked_C(cache):
    """C chosen on the complete labelled set by repeated stratified 5-fold
    AUC (5 repeats), as in the exploratory ensemble's tuning step."""
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=SEED)
    scores = np.zeros(len(C_GRID))
    for tr, te in rskf.split(X, y):
        Xa, Xb = cache.get(tr, te)
        sel = SelectKBest(f_classif, k=9).fit(Xa, y[tr])
        Xa, Xb = sel.transform(Xa), sel.transform(Xb)
        for j, C in enumerate(C_GRID):
            m = LogisticRegression(class_weight="balanced", max_iter=2000, C=C).fit(Xa, y[tr])
            scores[j] += roc_auc_score(y[te], m.predict_proba(Xb)[:, 1])
    return float(C_GRID[int(np.argmax(scores))])


def run_lr_config(split, leak_scale, leak_select, leak_tune, C_full, cache):
    Xs = cache.full() if leak_scale else X
    mask = None
    if leak_select:
        mask = SelectKBest(f_classif, k=9).fit(cache.full(), y).get_support()
    folds = list(LeaveOneOut().split(X)) if split == "variant" else \
        list(LeaveOneGroupOut().split(X, y, g))
    preds = np.zeros(len(y))
    chosen = []
    for tr, te in folds:
        Xsrc = Xs[:, mask] if mask is not None else Xs
        if leak_tune:
            C = C_full
        else:
            # inner preprocessing must also be fitted inside the inner folds
            # when scaling is in fold; with a leaked selection mask the inner
            # data are restricted to the selected columns of the cached transform
            if leak_scale:
                sc = inner_auc_all(tr, Xsrc, cache, True, leak_select, split == "gene")
            else:
                sc = inner_auc_all(tr, None, cache if mask is None else MaskedCache(cache, mask),
                                   False, leak_select, split == "gene")
            C = float(C_GRID[int(np.argmax(sc))])
        chosen.append(C)
        if leak_scale:
            Xa, Xb = Xsrc[tr], Xsrc[te]
        else:
            Xa, Xb = cache.get(tr, te)
            if mask is not None:
                Xa, Xb = Xa[:, mask], Xb[:, mask]
        if not leak_select and Xa.shape[1] > 9:
            sel = SelectKBest(f_classif, k=9).fit(Xa, y[tr])
            Xa, Xb = sel.transform(Xa), sel.transform(Xb)
        m = LogisticRegression(class_weight="balanced", max_iter=2000, C=C).fit(Xa, y[tr])
        preds[te] = m.predict_proba(Xb)[:, 1]
    d = auc_decompose(y, preds, g)
    return dict(split=split,
                scaling="full set" if leak_scale else "in fold",
                selection="full set" if leak_select else "in fold",
                tuning="full set" if leak_tune else "in fold",
                auc=round(roc_auc_score(y, preds), 4),
                ap=round(average_precision_score(y, preds), 4),
                auc_within=round(d["auc_within"], 4), auc_between=round(d["auc_between"], 4),
                median_C=float(np.median(chosen))), preds


class MaskedCache:
    """View of a FoldCache restricted to a fixed column mask."""

    def __init__(self, cache, mask):
        self.cache, self.mask = cache, mask

    def get(self, a, b):
        Xa, Xb = self.cache.get(a, b)
        return Xa[:, self.mask], Xb[:, self.mask]


if __name__ == "__main__":
    # ── logistic regression, 2^4 ─────────────────────────────────────────
    cache = FoldCache(X)
    C_full = leaked_C(cache)
    print(f"C selected on the complete labelled set: {C_full:.4g}", flush=True)
    grid = list(itertools.product(["variant", "gene"], [True, False], [True, False], [True, False]))
    res = [run_lr_config(s, a, b, c, C_full, cache) for s, a, b, c in grid]
    rows, predcols = [], {}
    for i, (r, p) in enumerate(res):
        r = dict(config=f"L{i + 1}", **r)
        rows.append(r)
        predcols[r["config"]] = p
        print(r, flush=True)
    lr16 = pd.DataFrame(rows)
    lr16.to_csv(out("ablation_lr16.csv"), index=False)
    pd.DataFrame(predcols).assign(Gene_Name=g, label=y).to_csv(
        out("ablation_predictions.csv"), index=False)

    # factorial effects (Yates contrasts) on AUC: +1 = stage sees held-out data
    F = pd.DataFrame({
        "split": np.where(lr16.split == "variant", 1, -1),
        "scaling": np.where(lr16.scaling == "full set", 1, -1),
        "selection": np.where(lr16.selection == "full set", 1, -1),
        "tuning": np.where(lr16.tuning == "full set", 1, -1)})
    eff = []
    for k in range(1, 5):
        for combo in itertools.combinations(F.columns, k):
            sign = np.prod(F[list(combo)].values, axis=1)
            eff.append(dict(effect=" x ".join(combo), order=k,
                            auc_effect=round(float(np.mean(sign * lr16.auc) * 2), 4)))
    eff = pd.DataFrame(eff)
    eff.to_csv(out("ablation_effects.csv"), index=False)
    print(eff.sort_values("auc_effect", key=abs, ascending=False).head(8).to_string(index=False))

    # ── three model families, 2^3 ────────────────────────────────────────
    def one(split, ls, lsel, model):
        p, a, ap = cross_validate(X, y, g, split, model, k=9, leak_scale=ls, leak_select=lsel,
                                  cache=cache)
        d = auc_decompose(y, p, g)
        return (split, ls, lsel, model, round(a, 4), round(ap, 4),
                round(d["auc_within"], 4), round(d["auc_between"], 4))

    jobs = [(s, a, b, m) for s in ["variant", "gene"] for a in [True, False]
            for b in [True, False] for m in ["lr", "rf", "gb"]]
    res3 = [one(*j) for j in jobs]
    rows = {}
    for s, a, b, m, auc, ap, aw, ab in res3:
        key = (s, a, b)
        rows.setdefault(key, dict(split=s, scaling="full set" if a else "in fold",
                                  selection="full set" if b else "in fold"))
        rows[key].update({m: auc, f"{m}_ap": ap, f"{m}_within": aw, f"{m}_between": ab})
    full = pd.DataFrame(list(rows.values()))
    full.insert(0, "config", [f"C{i + 1}" for i in range(len(full))])
    full.to_csv(out("ablation_full.csv"), index=False)
    print(full.to_string(index=False))
