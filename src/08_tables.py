#!/usr/bin/env python3
"""
Step 8 - supplementary tables that combine several analyses.

  TableS1_discrimination.csv   AUC, AP, gene-cluster and naive bootstrap
                               intervals, width ratio, Obuchowski SE and the
                               within/between decomposition, for raw scores
                               and gene-held-out models
  TableS2_paired_tests.csv     DeLong versus clustered (Obuchowski) paired tests
"""
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

from depbench.common import (AA_COLS, ALL_COLS, DEP_COLS, N_BOOT, SEED,
                             cross_validate, gene_prevalence, out)
from depbench.metrics import (auc_decompose, cluster_bootstrap, delong_test,
                              delong_var_single, fast_auc, obuchowski_test,
                              obuchowski_var, percentile_ci)

lab = pd.read_csv(out("cohort_full.csv"))
y, g = lab.label.values, lab.Gene_Name.values
ens = pd.read_csv(out("ensemble_predictions.csv"))
assert (ens.Protein_Change.values == lab.Protein_Change.values).all()

scores = {
    "AlphaMissense (variant score)": lab.AM_Specific_Score.values,
    "AlphaMissense (residue mean)": lab.AM_Mean_Residue_Score.values,
    "DCS": lab.DCS_Score.values, "RCS": lab.RCS_Score.values,
    "Normalised DCS": lab.Normalized_DCS.values,
    "Conservation (reconstructed)": lab.conservation.values,
    "Normalised mobility (reconstructed)": lab.Norm_Dynamics.values,
    "|delta hydrophobicity|": lab.absHydro.values, "|delta volume|": lab.absSize.values,
    "|delta charge|": lab.absCharge.values,
    "Gene prevalence (variant excluded)": gene_prevalence(lab, True),
    "Exploratory ensemble (variant-level, as reported)": ens.ensemble_reported.values,
}
fitted = {
    "Full 19-feature model (C8)": (ALL_COLS, 9),
    "Position-level scores only": (DEP_COLS, None),
    "Substitution descriptors only": (AA_COLS, 9),
    "AlphaMissense only": (["AM_Specific_Score"], None),
    "AlphaMissense + position-level scores": (["AM_Specific_Score"] + DEP_COLS, None),
    "AlphaMissense + full feature set": (["AM_Specific_Score"] + ALL_COLS, 9),
}
for name, (cols, k) in fitted.items():
    scores[name + " †"] = cross_validate(lab[cols].values, y, g, "gene", "lr", k=k)[0]

rows = []
for name, s in scores.items():
    s = np.asarray(s, dtype=float)
    stat = lambda yb, idx, gb, s=s: fast_auc(yb, s[idx])
    cb = cluster_bootstrap(y, g, stat, B=N_BOOT, seed=SEED, cluster=True)
    nb = cluster_bootstrap(y, g, stat, B=N_BOOT, seed=SEED, cluster=False)
    clo, chi = percentile_ci(cb)
    nlo, nhi = percentile_ci(nb)
    d = auc_decompose(y, s, g)
    _, C = obuchowski_var(y, [s], g)
    _, vdl = delong_var_single(y, s)
    rows.append(dict(score=name, auc=round(d["auc"], 4), ap=round(average_precision_score(y, s), 4),
                     auc_within=round(d["auc_within"], 4), auc_between=round(d["auc_between"], 4),
                     cluster_ci=f"{clo:.3f}-{chi:.3f}", cluster_se=round(cb.std(ddof=1), 3),
                     naive_ci=f"{nlo:.3f}-{nhi:.3f}", naive_se=round(nb.std(ddof=1), 3),
                     width_ratio=round((chi - clo) / (nhi - nlo), 2),
                     obuchowski_se=round(float(np.sqrt(C[0, 0])), 3),
                     delong_se=round(float(np.sqrt(vdl)), 3)))
t1 = pd.DataFrame(rows)
t1.to_csv(out("TableS1_discrimination.csv"), index=False)
print(t1.to_string(index=False))

pairs = [
    ("AlphaMissense (variant score)", "Exploratory ensemble (variant-level, as reported)"),
    ("Exploratory ensemble (variant-level, as reported)", "Gene prevalence (variant excluded)"),
    ("AlphaMissense (variant score)", "Full 19-feature model (C8) †"),
    ("AlphaMissense (variant score)", "AlphaMissense only †"),
    ("AlphaMissense only †", "AlphaMissense + position-level scores †"),
    ("AlphaMissense only †", "AlphaMissense + full feature set †"),
    ("Full 19-feature model (C8) †", "Substitution descriptors only †"),
    ("Full 19-feature model (C8) †", "Position-level scores only †"),
    ("Substitution descriptors only †", "Position-level scores only †"),
    ("DCS", "Conservation (reconstructed)"),
    ("DCS", "RCS"),
    ("AlphaMissense (residue mean)", "DCS"),
    ("AlphaMissense (variant score)", "AlphaMissense (residue mean)"),
]
rows = []
for a, b in pairs:
    a1, a2, zd, pdl = delong_test(y, scores[a], scores[b])
    _, _, zo, po, se = obuchowski_test(y, scores[a], scores[b], g)
    rows.append(dict(comparison=f"{a} vs {b}", auc_1=round(a1, 3), auc_2=round(a2, 3),
                     difference=round(a1 - a2, 3), delong_z=round(zd, 2), delong_p=float(f"{pdl:.3g}"),
                     clustered_z=round(zo, 2), clustered_p=float(f"{po:.3g}"),
                     clustered_se_diff=round(se, 3)))
t2 = pd.DataFrame(rows)
t2.to_csv(out("TableS2_paired_tests.csv"), index=False)
print(t2.to_string(index=False))
