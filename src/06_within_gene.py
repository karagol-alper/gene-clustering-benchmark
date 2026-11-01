#!/usr/bin/env python3
"""
Step 6 - within-gene discrimination of each raw score.

For every score: pooled within-gene AUC (covariate-adjusted AUC with gene as
stratum), per-gene AUCs, an exact stratified permutation test, a stratified
bootstrap interval (variants resampled within gene and class), the between-
gene AUC, and conditional logistic regression stratified by gene (odds ratio
per cohort standard deviation; and per within-protein percentile unit as a
rank-based check).
"""
import warnings

import numpy as np
import pandas as pd
from statsmodels.discrete.conditional_models import ConditionalLogit

from depbench.common import N_BOOT, SEED, gene_prevalence, out
from depbench.metrics import (auc_decompose, cluster_bootstrap, fast_auc,
                              per_gene_auc, percentile_ci, within_gene_perm_p)

warnings.filterwarnings("ignore")
lab = pd.read_csv(out("cohort_full.csv"))
y, g = lab.label.values, lab.Gene_Name.values
SCORES = {
    "AlphaMissense (variant)": "AM_Specific_Score",
    "AlphaMissense (residue mean)": "AM_Mean_Residue_Score",
    "DCS": "DCS_Score", "RCS": "RCS_Score", "Normalised DCS": "Normalized_DCS",
    "Raw DCS (before rescaling)": "Raw_DCS", "Raw RCS (before rescaling)": "Raw_RCS",
    "Conservation": "conservation", "Normalised mobility": "Norm_Dynamics",
    "Mean RMSF": "Avg_RMSF",
    "|delta hydrophobicity|": "absHydro", "|delta volume|": "absSize", "|delta charge|": "absCharge",
}
mixed = [x for x in np.unique(g) if 0 < y[g == x].sum() < (g == x).sum()]
rng = np.random.default_rng(SEED)
strata = [np.where((g == x) & (y == c))[0] for x in mixed for c in (0, 1)]


def strat_boot(s, B=N_BOOT):
    out_ = []
    for _ in range(B):
        idx = np.concatenate([rng.choice(ix, len(ix)) for ix in strata])
        out_.append(auc_decompose(y[idx], s[idx], g[idx])["auc_within"])
    return percentile_ci(np.array(out_))


rows = []
for name, col in SCORES.items():
    s = lab[col].values.astype(float)
    d = auc_decompose(y, s, g)
    pg = per_gene_auc(y, s, g)
    ow, p1, p2, _ = within_gene_perm_p(y, s, g, n_perm=20000, seed=SEED)
    lo, hi = strat_boot(s)
    clo, chi = percentile_ci(cluster_bootstrap(y, g, lambda yb, idx, gb, s=s: fast_auc(yb, s[idx]),
                                               B=N_BOOT, seed=SEED))
    z = (s - s.mean()) / s.std()
    cl = ConditionalLogit(y, z[:, None], groups=g).fit(disp=0)
    r = dict(score=name, auc=round(d["auc"], 4), auc_cluster_ci=f"{clo:.3f}-{chi:.3f}",
             auc_between=round(d["auc_between"], 4), auc_within=round(d["auc_within"], 4),
             within_ci=f"{lo:.3f}-{hi:.3f}", within_p_two_sided=round(p2, 4),
             within_p_one_sided=round(p1, 4),
             **{f"auc_{k}": round(v, 3) for k, v in pg.items()},
             clr_or_per_sd=round(float(np.exp(cl.params[0])), 3),
             clr_or_ci=f"{np.exp(cl.conf_int()[0, 0]):.2f}-{np.exp(cl.conf_int()[0, 1]):.2f}",
             clr_p=round(float(cl.pvalues[0]), 4))
    pc = f"{col}__pct_ref"
    if pc in lab:
        cl2 = ConditionalLogit(y, lab[pc].values[:, None] * 10, groups=g).fit(disp=0)
        r.update(clr_rank_or_per_decile=round(float(np.exp(cl2.params[0])), 3),
                 clr_rank_p=round(float(cl2.pvalues[0]), 4))
    rows.append(r)
    print(r, flush=True)
tab = pd.DataFrame(rows)
tab.to_csv(out("within_gene_discrimination.csv"), index=False)

# per-gene composition table (Table S3)
comp = lab.groupby("Gene_Name").agg(n=("label", "size"), n_pathogenic=("label", "sum"),
                                    uniprot=("UniProt_ID", "first"),
                                    dcs_multiplier=("dcs_scale_factor", "first"),
                                    median_dcs=("DCS_Score", "median"),
                                    median_conservation=("conservation", "median"),
                                    median_am=("AM_Specific_Score", "median")).reset_index()
comp["n_benign"] = comp.n - comp.n_pathogenic
comp["pathogenic_fraction"] = (comp.n_pathogenic / comp.n).round(3)
comp["within_pairs"] = comp.n_pathogenic * comp.n_benign
for name, col in [("AM", "AM_Specific_Score"), ("DCS", "DCS_Score"), ("RCS", "RCS_Score"),
                  ("conservation", "conservation")]:
    pg = per_gene_auc(y, lab[col].values, g)
    comp[f"auc_{name}"] = comp.Gene_Name.map(pg).round(3)
comp = comp.sort_values("n", ascending=False)
comp.to_csv(out("TableS_per_gene.csv"), index=False)
print(comp.to_string(index=False))
