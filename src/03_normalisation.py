#!/usr/bin/env python3
"""
Step 3 - provenance, protein identity and cross-protein transfer of the
position-level scores.

  A  per-protein rescaling parameters implied by the published scores
  B  gene-identity recovery from each encoding (labelled cohort and the
     label-free set of all ClinVar-observed residues), residue-grouped CV,
     balanced accuracy and a label-permutation null
  C  gene-held-out transfer of each encoding, with the within/between
     decomposition of the pooled AUC
  D  incremental value over AlphaMissense under gene-held-out validation,
     with DeLong and clustered (Obuchowski) paired tests and a paired
     gene-cluster bootstrap of the AUC difference
"""
import json

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.model_selection import StratifiedGroupKFold

from depbench.common import (AA_COLS, ALL_COLS, DEP_COLS, N_BOOT, N_JOBS, SEED,
                             cross_validate, out)
from depbench.metrics import (auc_decompose, cluster_bootstrap, delong_test,
                              fast_auc, obuchowski_test, percentile_ci)

lab = pd.read_csv(out("cohort_full.csv"))
allv = pd.read_csv(out("clinvar_all_components.csv"))
maps = pd.read_csv(out("protein_scaling_maps.csv"))
y, g = lab.label.values, lab.Gene_Name.values

# transplant of the published per-protein DCS multipliers onto AlphaMissense
for df in (lab, allv):
    fac = df.dcs_scale_factor / lab.dcs_scale_factor.mean()
    df["AM_Mean_Residue_Score__xdcs"] = df.AM_Mean_Residue_Score * fac
    df["AM_Specific_Score__xdcs"] = df.AM_Specific_Score * fac

pct = lambda cols, t="ref": [f"{c}__pct_{t}" for c in cols]
ENC = {   # name: (columns, description)
    "published":       (DEP_COLS, "DCS, RCS, normalised DCS as published (per-protein min-max, 0-10)"),
    "unrescaled":      (["Raw_DCS", "Raw_RCS"], "conservation x mobility products before per-protein rescaling"),
    "components":      (["conservation", "Norm_Dynamics"], "conservation and normalised mobility"),
    "conservation":    (["conservation"], "conservation alone"),
    "mobility":        (["Norm_Dynamics"], "normalised mobility alone"),
    "rmsf":            (["Avg_RMSF"], "mean backbone RMSF (absolute, Angstrom)"),
    "pct_protein":     (pct(DEP_COLS), "within-protein percentile, label-free reference"),
    "pct_cohort":      (pct(DEP_COLS, "coh"), "within-gene percentile among labelled variants"),
    "z_cohort":        ([f"{c}__z_coh" for c in DEP_COLS], "within-gene z-score among labelled variants"),
    "ctr_cohort":      ([f"{c}__ctr_coh" for c in DEP_COLS], "within-gene centring among labelled variants"),
    "am_res":          (["AM_Mean_Residue_Score"], "AlphaMissense residue mean (absolute scale)"),
    "am_res_minmax":   (["AM_Mean_Residue_Score__mm_ref"], "AlphaMissense residue mean, per-protein min-max"),
    "am_res_xdcs":     (["AM_Mean_Residue_Score__xdcs"], "AlphaMissense residue mean x published DCS multiplier"),
    "am_res_pct":      (["AM_Mean_Residue_Score__pct_ref"], "AlphaMissense residue mean, within-protein percentile"),
    "substitution":    (AA_COLS, "16 substitution descriptors"),
}
results = {}

# ── A. provenance of the per-protein rescaling ──────────────────────────
coh_genes = lab.groupby("Protein_ID").Gene_Name.first()
nd_max = allv.groupby("Protein_ID").Norm_Dynamics.max()
prov = maps.assign(
    in_cohort=maps.Protein_ID.isin(coh_genes.index),
    rmsf_at_nd0=maps.rmsf_icpt, rmsf_at_nd1=maps.rmsf_slope + maps.rmsf_icpt,
    max_nd_clinvar_residues=maps.Protein_ID.map(nd_max))
prov.to_csv(out("TableS_provenance.csv"), index=False)
pc = prov[prov.in_cohort]
results["provenance"] = dict(
    dcs_mult_min=float(pc.dcs_slope.min()), dcs_mult_max=float(pc.dcs_slope.max()),
    rcs_mult_min=float(pc.rcs_slope.min()), rcs_mult_max=float(pc.rcs_slope.max()),
    rmsf_max_min=float(pc.rmsf_at_nd1.min()), rmsf_max_max=float(pc.rmsf_at_nd1.max()),
    n_cohort_proteins=int(len(pc)),
    n_max_nd_below_0_5=int((pc.max_nd_clinvar_residues < 0.5).sum()),
    median_max_nd=float(pc.max_nd_clinvar_residues.median()))
print("provenance:", results["provenance"])


# ── B. gene-identity recovery ───────────────────────────────────────────
def gene_id(df, cols, gene_col, pos_key, n_rep=5, n_perm=50, seed=SEED):
    Xg = df[cols].values
    yg = pd.factorize(df[gene_col])[0]
    grp = df[pos_key].values

    def score(labels, rep):
        cv = StratifiedGroupKFold(5, shuffle=True, random_state=seed + rep)
        pred = np.empty(len(labels), dtype=int)
        for tr, te in cv.split(Xg, labels, grp):
            rf = RandomForestClassifier(n_estimators=150, random_state=seed + rep, n_jobs=1)
            pred[te] = rf.fit(Xg[tr], labels[tr]).predict(Xg[te])
        return balanced_accuracy_score(labels, pred), accuracy_score(labels, pred)

    obs = np.mean([score(yg, r) for r in range(n_rep)], axis=0)
    rng = np.random.default_rng(seed)
    null = np.array([score(rng.permutation(yg), 100 + b)[0] for b in range(n_perm)])
    return dict(bacc=obs[0], acc=obs[1], null_mean=null.mean(), null_sd=null.std(ddof=1),
                p=(np.sum(null >= obs[0]) + 1) / (n_perm + 1), k=len(np.unique(yg)),
                majority=pd.Series(yg).value_counts(normalize=True).max())


big = [x for x, n in lab.Gene_Name.value_counts().items() if n >= 3]
sub = lab[lab.Gene_Name.isin(big)].copy()
sub["pos_key"] = sub.Protein_ID + "_" + sub.Residue_Pos.astype(str)
ref = allv.drop_duplicates(["Protein_ID", "Residue_Pos"]).copy()
ref["pos_key"] = ref.Protein_ID + "_" + ref.Residue_Pos.astype(str)
big_ref = [p for p, n in ref.Protein_ID.value_counts().items() if n >= 10]
ref = ref[ref.Protein_ID.isin(big_ref)]

gid_enc = ["published", "unrescaled", "components", "conservation", "mobility", "rmsf",
           "pct_protein", "am_res", "am_res_minmax", "am_res_xdcs", "am_res_pct", "substitution"]
jobs = [("labelled", e) for e in gid_enc] + \
       [("all_clinvar", e) for e in gid_enc if e not in ("substitution",)]


def run_gid(setname, e):
    cols = ENC[e][0]
    if setname == "labelled":
        return setname, e, gene_id(sub, cols, "Gene_Name", "pos_key")
    if any(c.endswith("__pct_coh") for c in cols):
        return setname, e, None
    return setname, e, gene_id(ref, cols, "Gene_Name", "pos_key")


rows = []
for setname, e, r in Parallel(n_jobs=N_JOBS)(delayed(run_gid)(s, e) for s, e in jobs):
    if r is None:
        continue
    rows.append(dict(set=setname, encoding=e, description=ENC[e][1],
                     **{k: round(float(v), 4) for k, v in r.items()}))
    print(rows[-1], flush=True)
gid = pd.DataFrame(rows)
gid.to_csv(out("gene_identity_recovery.csv"), index=False)
results["gene_identity"] = dict(n_genes_labelled=len(big), n_variants_labelled=int(len(sub)),
                                n_proteins_all=len(big_ref), n_residues_all=int(len(ref)))

# ── C. gene-held-out transfer of each encoding ──────────────────────────
tr_enc = ["published", "unrescaled", "components", "conservation", "mobility", "rmsf",
          "pct_protein", "pct_cohort", "z_cohort", "ctr_cohort",
          "am_res", "am_res_minmax", "am_res_xdcs", "am_res_pct", "substitution"]


def transfer(e):
    cols = ENC[e][0]
    k = 9 if len(cols) > 9 else None
    pg, a, ap = cross_validate(lab[cols].values, y, g, "gene", "lr", k=k)
    pv, av, apv = cross_validate(lab[cols].values, y, g, "variant", "lr", k=k)
    dg = auc_decompose(y, pg, g)
    return dict(encoding=e, description=ENC[e][1], n_features=len(cols),
                auc_gene=round(a, 4), ap_gene=round(ap, 4),
                within_gene=round(dg["auc_within"], 4), between_gene=round(dg["auc_between"], 4),
                auc_variant=round(av, 4), ap_variant=round(apv, 4)), pg


out_c = Parallel(n_jobs=N_JOBS)(delayed(transfer)(e) for e in tr_enc)
trans = pd.DataFrame([r for r, _ in out_c])
trans.to_csv(out("transfer_by_encoding.csv"), index=False)
print(trans.to_string(index=False))
logo_preds = {r["encoding"]: p for r, p in out_c}

# ── D. incremental value over AlphaMissense (gene-held-out) ─────────────
AM = ["AM_Specific_Score"]
INC = {
    "AM alone": AM,
    "AM + position scores (published)": AM + DEP_COLS,
    "AM + position scores (unrescaled)": AM + ["Raw_DCS", "Raw_RCS"],
    "AM + conservation": AM + ["conservation"],
    "AM + position scores (within-protein percentile)": AM + pct(DEP_COLS),
    "AM + position scores (within-gene percentile, labelled)": AM + pct(DEP_COLS, "coh"),
    "AM + position scores (within-gene z, labelled)": AM + [f"{c}__z_coh" for c in DEP_COLS],
    "AM + position scores (within-gene centred, labelled)": AM + [f"{c}__ctr_coh" for c in DEP_COLS],
    "AM + substitution descriptors": AM + AA_COLS,
    "AM + full feature set": AM + ALL_COLS,
}


def inc(name):
    cols = INC[name]
    k = 9 if len(cols) > 9 else None
    p, a, ap = cross_validate(lab[cols].values, y, g, "gene", "lr", k=k)
    return name, p, a, ap


inc_res = Parallel(n_jobs=N_JOBS)(delayed(inc)(n) for n in INC)
preds = {n: p for n, p, _, _ in inc_res}
base = preds["AM alone"]
rows = []
for n, p, a, ap in inc_res:
    d = auc_decompose(y, p, g)
    r = dict(model=n, n_features=len(INC[n]), feature_selection="k=9" if len(INC[n]) > 9 else "none",
             auc=round(a, 4), ap=round(ap, 4), within=round(d["auc_within"], 4),
             between=round(d["auc_between"], 4))
    if n != "AM alone":
        _, _, zd, pd_ = delong_test(y, p, base)
        _, _, zo, po, se = obuchowski_test(y, p, base, g)

        def diff(yb, idx, gb, p=p):
            return fast_auc(yb, p[idx]) - fast_auc(yb, base[idx])
        boot = cluster_bootstrap(y, g, diff, B=N_BOOT, seed=SEED)
        lo, hi = percentile_ci(boot)
        r.update(delta=round(a - fast_auc(y, base), 4), delong_z=round(zd, 3), delong_p=round(pd_, 4),
                 clustered_z=round(zo, 3), clustered_p=round(po, 4),
                 delta_cluster_ci_lo=round(float(lo), 4), delta_cluster_ci_hi=round(float(hi), 4))
    rows.append(r)
incr = pd.DataFrame(rows)
incr.to_csv(out("incremental_over_am.csv"), index=False)
print(incr.to_string(index=False))
pd.DataFrame({**{f"logo::{k}": v for k, v in logo_preds.items()},
              **{f"inc::{k}": v for k, v in preds.items()}}).assign(
    Gene_Name=g, label=y).to_csv(out("normalisation_predictions.csv"), index=False)
json.dump(results, open(out("normalisation.json"), "w"), indent=1, default=float)
