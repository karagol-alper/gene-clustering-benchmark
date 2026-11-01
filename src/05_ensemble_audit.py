#!/usr/bin/env python3
"""
Step 5 - audit of the tuned five-model ensemble (AUC 0.756).

The ensemble's out-of-fold predictions (leave-one-variant-out) for its five
base learners are distributed in loo_predictions.csv.  Its blending weights
were chosen by 300 trials of Bayesian optimisation maximising the
leave-one-variant-out AUC itself, so the reported value contains a final
in-sample selection step.  This script:
  * recovers the blending weights from the distributed predictions
  * re-estimates the blend with the weight search nested inside the
    leave-one-out loop (weights chosen without the held-out variant)
  * compares the ensemble with the gene-prevalence baseline (DeLong and
    clustered tests) and decomposes both into within/between-gene parts
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

from depbench.common import DATA_DIR, N_BOOT, SEED, gene_prevalence, out
from depbench.metrics import (auc_decompose, cluster_bootstrap, delong_test,
                              fast_auc, obuchowski_test, percentile_ci,
                              within_gene_perm_p)

lab = pd.read_csv(out("cohort_full.csv"))
loo = pd.read_csv(f"{DATA_DIR}/loo_predictions.csv")
m = lab.merge(loo.drop(columns=["Protein_ID", "ClinVar_Significance", "DCS_Score", "RCS_Score"]),
              on=["Gene_Name", "Protein_Change"], how="left", validate="1:1")
assert m.pred_ensemble.notna().all() and (m.true_label == m.label).all()
y, g = m.label.values, m.Gene_Name.values
BASES = ["pred_XGBoost", "pred_LightGBM", "pred_CatBoost", "pred_RandomForest", "pred_Logistic"]
P = m[BASES].values
ens = m.pred_ensemble.values
res = {}

# recovered weights (the ensemble is an exact convex combination)
w, *_ = np.linalg.lstsq(P, ens, rcond=None)
res["recovered_weights"] = {b.replace("pred_", ""): round(float(x), 4) for b, x in zip(BASES, w)}
res["recovery_max_abs_error"] = float(np.abs(P @ w - ens).max())
res["auc_reported"] = fast_auc(y, ens)
res["auc_base_models"] = {b.replace("pred_", ""): round(fast_auc(y, P[:, j]), 4) for j, b in enumerate(BASES)}
res["auc_equal_weights"] = fast_auc(y, P.mean(1))

# candidate weight vectors: vertices, equal weights and Dirichlet(1) draws
rng = np.random.default_rng(SEED)
W = np.vstack([np.eye(5), np.full(5, 0.2), rng.dirichlet(np.ones(5), size=3000)])
B = P @ W.T                                            # n x n_candidates
auc_all = np.array([fast_auc(y, B[:, j]) for j in range(B.shape[1])])
res["auc_in_sample_best"] = float(auc_all.max())
res["in_sample_best_weights"] = W[auc_all.argmax()].round(3).tolist()

# nested: weights chosen on the other 87 variants
nested = np.empty(len(y))
for i in range(len(y)):
    keep = np.arange(len(y)) != i
    a = np.array([fast_auc(y[keep], B[keep, j]) for j in range(B.shape[1])])
    nested[i] = B[i, a.argmax()]
res["auc_nested_blend"] = fast_auc(y, nested)

# nested choice of the single best base learner
nested_sel = np.empty(len(y))
for i in range(len(y)):
    keep = np.arange(len(y)) != i
    a = [fast_auc(y[keep], P[keep, j]) for j in range(5)]
    nested_sel[i] = P[i, int(np.argmax(a))]
res["auc_nested_best_single"] = fast_auc(y, nested_sel)

# comparison with the gene-prevalence lookup
prev = gene_prevalence(m, leave_one_out=True)
for name, s in [("ensemble", ens), ("nested_blend", nested), ("prevalence_loo", prev)]:
    d = auc_decompose(y, s, g)
    lo, hi = percentile_ci(cluster_bootstrap(y, g, lambda yb, idx, gb, s=s: fast_auc(yb, s[idx]),
                                             B=N_BOOT, seed=SEED))
    res[f"{name}_decomposition"] = {k: float(v) for k, v in d.items()}
    res[f"{name}_cluster_ci"] = [float(lo), float(hi)]
_, _, zd, pd_ = delong_test(y, ens, prev)
_, _, zo, po, se = obuchowski_test(y, ens, prev, g)
res["ens_vs_prev"] = dict(delta=float(fast_auc(y, ens) - fast_auc(y, prev)), delong_z=float(zd),
                          delong_p=float(pd_), clustered_z=float(zo), clustered_p=float(po),
                          clustered_se=float(se))
res["spearman_ens_prevalence_in_sample"] = float(stats.spearmanr(ens, gene_prevalence(m)).statistic)
res["spearman_ens_prevalence_loo"] = float(stats.spearmanr(ens, prev).statistic)
ow, p1, p2, nm = within_gene_perm_p(y, ens, g, n_perm=10000, seed=SEED)
res["ensemble_within_gene"] = dict(auc_within=float(ow), p_one_sided=float(p1), p_two_sided=float(p2))
json.dump(res, open(out("ensemble_audit.json"), "w"), indent=1, default=float)
pd.DataFrame({"Gene_Name": g, "Protein_Change": m.Protein_Change, "label": y,
              "ensemble_reported": ens, "ensemble_nested_blend": nested,
              "gene_prevalence_loo": prev}).to_csv(out("ensemble_predictions.csv"), index=False)
for k, v in res.items():
    print(k, v)
