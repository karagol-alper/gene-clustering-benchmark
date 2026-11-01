#!/usr/bin/env python3
"""
Step 3b - sensitivity of within-gene percentile encodings to the ranking
convention.  Mid-rank percentiles, (rank - 1/2)/n, have mean 1/2 in every gene;
rank/n (pandas rank(pct=True)) has mean (n + 1)/2n and so encodes gene size.
"""
import json

import numpy as np
import pandas as pd

from depbench.common import DEP_COLS, cross_validate, out
from depbench.metrics import auc_decompose, fast_auc

lab = pd.read_csv(out("cohort_full.csv"))
y, g = lab.label.values, lab.Gene_Name.values
am = lab.AM_Specific_Score.values
res = {}
for name, P in [("midrank", lab[[f"{c}__pct_coh" for c in DEP_COLS]].values),
                ("rank_over_n", lab.groupby("Gene_Name")[DEP_COLS].rank(pct=True).values)]:
    p_am, a_am, _ = cross_validate(np.c_[am, P], y, g, "gene", "lr", k=None)
    p_alone, a_alone, _ = cross_validate(P, y, g, "gene", "lr", k=None)
    d = auc_decompose(y, p_am, g)
    res[name] = dict(am_plus_pct=a_am, pct_alone=a_alone, within=d["auc_within"],
                     between=d["auc_between"],
                     mean_by_size={int(k): float(v) for k, v in
                                   pd.Series(P[:, 0]).groupby(lab.n_in_gene.values).mean().items()})
base = cross_validate(am[:, None], y, g, "gene", "lr", k=None)[1]
res["am_alone"] = base
for k in ("midrank", "rank_over_n"):
    res[k]["delta_vs_am"] = res[k]["am_plus_pct"] - base
json.dump(res, open(out("rank_convention.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
