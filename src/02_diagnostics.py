#!/usr/bin/env python3
"""
Step 2 - cluster diagnostics.

  * gene-class association: chi-square with a Monte-Carlo exact p-value,
    Cramer's V, and the between-gene SD of log-odds from a random-intercept
    logistic model (profile-likelihood interval)
  * composition of positive-negative pairs (within versus between genes)
  * gene-prevalence baselines: in sample, leave-one-variant-out, and
    leave-one-gene-out (the prior-shift floor of pooled gene-held-out AUC)
  * permutation nulls (global and within-gene label permutation) for the
    variant-level pipeline (C1) and the gene-level pipeline (C8), and for
    the position-level scores evaluated gene-held-out, with a third null in
    which each gene's position-level scores are replaced by residues drawn
    from a randomly assigned donor protein other than its own (gene
    fingerprint, no signal).  Upper-tail, lower-tail and two-sided
    (twice the smaller tail) Monte Carlo p-values are recorded.
"""
import json

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats

from depbench.common import (ALL_COLS, DEP_COLS, N_JOBS, N_PERM, SEED,
                             FoldCache, cross_validate, gene_prevalence,
                             gene_prevalence_logo, out)
from depbench.glmm import fit_ri_logit
from depbench.metrics import auc_decompose, fast_auc, mc_se, perm_p

lab = pd.read_csv(out("cohort_full.csv"))
allv = pd.read_csv(out("clinvar_all_components.csv"))
y, g = lab.label.values, lab.Gene_Name.values
X = lab[ALL_COLS].values.astype(float)
genes = np.unique(g)
res = {}

# ── gene-class association ────────────────────────────────────────────────
tab = pd.crosstab(lab.Gene_Name, lab.label)
chi2, p_asym, dof, expected = stats.chi2_contingency(tab)
rng = np.random.default_rng(SEED)
null = np.empty(100_000)
# Monte-Carlo exact test: permute labels, recompute chi-square from counts
gi = pd.factorize(g)[0]
n_g = np.bincount(gi)
exp_p = n_g * y.mean()
exp_n = n_g * (1 - y.mean())
for b in range(len(null)):
    yp = rng.permutation(y)
    o_p = np.bincount(gi, weights=yp, minlength=len(n_g))
    o_n = n_g - o_p
    null[b] = np.sum((o_p - exp_p) ** 2 / exp_p + (o_n - exp_n) ** 2 / exp_n)
p_mc = (np.sum(null >= chi2 - 1e-9) + 1) / (len(null) + 1)
cramer_v = np.sqrt(chi2 / len(y))
ri = fit_ri_logit(y, g)
res["association"] = dict(chi2=round(chi2, 2), df=int(dof), p_asymptotic=float(p_asym),
                          p_monte_carlo=float(p_mc), n_mc=len(null),
                          frac_expected_lt5=float((expected < 5).mean()),
                          cramers_v=round(cramer_v, 3), **{f"ri_{k}": v for k, v in ri.items()})
print("association:", res["association"])

# ── pair composition ────────────────────────────────────────────────────
d0 = auc_decompose(y, np.zeros(len(y)), g)
per_gene_pairs = {gn: int(((g == gn) & (y == 1)).sum() * ((g == gn) & (y == 0)).sum())
                  for gn in genes if 0 < y[g == gn].sum() < (g == gn).sum()}
res["pairs"] = dict(n_pairs=d0["n_pairs"], n_within=d0["n_within"], n_between=d0["n_between"],
                    w_within=d0["w_within"], per_gene=per_gene_pairs,
                    n_mixed_genes=len(per_gene_pairs),
                    n_single_class_genes=int(len(genes) - len(per_gene_pairs)),
                    n_variants_in_single_class_genes=int(sum((g == gn).sum() for gn in genes
                                                             if gn not in per_gene_pairs)))
print("pairs:", res["pairs"])

# ── gene-prevalence baselines ───────────────────────────────────────────
prev = {"in_sample": gene_prevalence(lab), "variant_loo": gene_prevalence(lab, True),
        "gene_loo": gene_prevalence_logo(lab)}
res["prevalence"] = {k: {kk: (float(vv) if not isinstance(vv, int) else vv)
                         for kk, vv in auc_decompose(y, v, g).items()} for k, v in prev.items()}
pd.DataFrame(prev).assign(Gene_Name=g, label=y).to_csv(out("gene_prevalence_scores.csv"), index=False)
for k, v in res["prevalence"].items():
    print(f"gene prevalence {k:12s} AUC {v['auc']:.3f} within {v['auc_within']:.3f} between {v['auc_between']:.3f}")


# ── permutation nulls ───────────────────────────────────────────────────
def within_perm(yy, rng_):
    yp = yy.copy()
    for gn in genes:
        idx = np.where(g == gn)[0]
        yp[idx] = rng_.permutation(yy[idx])
    return yp


ref = allv.drop_duplicates(["Protein_ID", "Residue_Pos"])
donor_ids = [p for p, n in ref.Protein_ID.value_counts().items() if n >= 5]
donor_tab = {p: ref.loc[ref.Protein_ID == p, DEP_COLS].values for p in donor_ids}


own_protein = lab.groupby("Gene_Name").Protein_ID.first().reindex(genes).values


def donor_features(rng_):
    """Each gene receives residues of a different protein (never its own)."""
    Xd = np.empty((len(y), len(DEP_COLS)))
    while True:
        pick = rng_.choice(donor_ids, size=len(genes), replace=False)
        if not np.any(pick == own_protein):
            break
    for gn, p in zip(genes, pick):
        idx = np.where(g == gn)[0]
        rows = donor_tab[p][rng_.integers(0, len(donor_tab[p]), len(idx))]
        Xd[idx] = rows
    return Xd


def run_null(kind, design, seed, n):
    """Null replicates.  Label-independent preprocessing is cached per split
    and the penalised logistic model is solved with the vectorised Newton
    solver (identical objective); observed values are recomputed the same way."""
    rng_ = np.random.default_rng(seed)
    cX, cP = FoldCache(X), FoldCache(lab[DEP_COLS].values)
    outv = []
    for _ in range(n):
        if design == "C1":
            yp = rng_.permutation(y) if kind == "global" else within_perm(y, rng_)
            p, _, _ = cross_validate(X, yp, g, "variant", "lr", k=9, leak_scale=True,
                                     leak_select=True, cache=cX, engine="fast")
        elif design == "C8":
            yp = rng_.permutation(y) if kind == "global" else within_perm(y, rng_)
            p, _, _ = cross_validate(X, yp, g, "gene", "lr", k=9, cache=cX, engine="fast")
        elif design == "POSV":  # position-level scores, variant-level
            yp = rng_.permutation(y) if kind == "global" else within_perm(y, rng_)
            p, _, _ = cross_validate(lab[DEP_COLS].values, yp, g, "variant", "lr", k=None,
                                     cache=cP, engine="fast")
        else:  # position-level scores, gene-held-out
            if kind == "donor":
                yp = y
                Xd = donor_features(rng_)
                p, _, _ = cross_validate(Xd, y, g, "gene", "lr", k=None, engine="fast")
            else:
                yp = rng_.permutation(y) if kind == "global" else within_perm(y, rng_)
                p, _, _ = cross_validate(lab[DEP_COLS].values, yp, g, "gene", "lr", k=None,
                                         cache=cP, engine="fast")
        d = auc_decompose(yp, p, g)
        outv.append((d["auc"], d["auc_within"], d["auc_between"]))
    return outv


observed = {
    "C1": cross_validate(X, y, g, "variant", "lr", k=9, leak_scale=True, leak_select=True,
                         engine="fast")[0],
    "C8": cross_validate(X, y, g, "gene", "lr", k=9, engine="fast")[0],
    "POS": cross_validate(lab[DEP_COLS].values, y, g, "gene", "lr", k=None, engine="fast")[0],
    "POSV": cross_validate(lab[DEP_COLS].values, y, g, "variant", "lr", k=None, engine="fast")[0],
}
obs = {k: auc_decompose(y, v, g) for k, v in observed.items()}
print({k: round(v["auc"], 4) for k, v in obs.items()})
plan = [("C1", "global", N_PERM), ("C1", "within_gene", N_PERM),
        ("C8", "global", N_PERM), ("C8", "within_gene", N_PERM),
        ("POS", "global", N_PERM), ("POS", "within_gene", N_PERM),
        ("POS", "donor", N_PERM // 2),
        ("POSV", "global", N_PERM), ("POSV", "within_gene", N_PERM)]
dists, summary = {}, []
for i, (design, kind, n) in enumerate(plan):
    chunks = [n // N_JOBS + (1 if j < n % N_JOBS else 0) for j in range(N_JOBS)]
    parts = Parallel(n_jobs=N_JOBS)(delayed(run_null)(kind, design, SEED * 100 + i * 10 + j, c)
                                    for j, c in enumerate(chunks))
    a = np.array([r for part in parts for r in part])
    dists[f"{design}_{kind}"] = a[:, 0]
    o = obs[design]["auc"]
    p_hi = perm_p(a[:, 0], o, greater=True)
    p_lo = perm_p(a[:, 0], o, greater=False)
    p_two = min(1.0, 2 * min(p_hi, p_lo))
    summary.append(dict(design=design, null=kind, n=len(a), observed=round(o, 4),
                        mean=round(a[:, 0].mean(), 4), sd=round(a[:, 0].std(ddof=1), 4),
                        pct2_5=round(np.percentile(a[:, 0], 2.5), 4),
                        pct97_5=round(np.percentile(a[:, 0], 97.5), 4),
                        p_upper=round(p_hi, 4), p_lower=round(p_lo, 4), p_two_sided=round(p_two, 4),
                        mc_se_upper=round(mc_se(p_hi, len(a)), 4),
                        mean_within=round(np.nanmean(a[:, 1]), 4),
                        mean_between=round(np.nanmean(a[:, 2]), 4)))
    print(summary[-1], flush=True)
pd.DataFrame({k: pd.Series(v) for k, v in dists.items()}).to_csv(
    out("permutation_distributions.csv"), index=False)
pd.DataFrame(summary).to_csv(out("permutation_summary.csv"), index=False)
res["observed"] = {k: {kk: float(vv) for kk, vv in v.items()} for k, v in obs.items()}
json.dump(res, open(out("diagnostics.json"), "w"), indent=1, default=float)
