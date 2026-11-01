#!/usr/bin/env python3
"""
Consistency audit: recompute key quantities directly from results/ (not via
the manuscript generator) and confirm that each appears in the manuscript
text; check for unresolved placeholders; check citation order and coverage.
Exits non-zero on any failure.
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd

from depbench.common import ROOT, out

MS = open(os.path.join(ROOT, "manuscript", "manuscript.md")).read()
SUPP = open(os.path.join(ROOT, "manuscript", "supplement.md")).read()
ok, bad = [], []


def expect(label, value, text=MS):
    good = value in text
    (ok if good else bad).append(f"{'PASS' if good else 'FAIL'}  {label}: '{value}'")


f3 = lambda x: f"{float(x):.3f}"
coh = pd.read_csv(out("cohort_full.csv"))
expect("cohort size", f"{len(coh)} ClinVar")
expect("genes", f"in {coh.Gene_Name.nunique()} genes")
expect("positives/negatives", f"{int(coh.label.sum())} P/LP, {int((coh.label == 0).sum())} B/LB")

d = json.load(open(out("diagnostics.json")))
expect("chi-square", f"χ² = {d['association']['chi2']:.1f}")
expect("between-gene pairs", f"{d['pairs']['n_between']:,} ({100 * (1 - d['pairs']['w_within']):.0f}%)")
expect("prevalence in sample", f3(d["prevalence"]["in_sample"]["auc"]))
expect("prevalence variant-LOO", f3(d["prevalence"]["variant_loo"]["auc"]))
expect("prevalence gene-LOO", f3(d["prevalence"]["gene_loo"]["auc"]))

ps = pd.read_csv(out("permutation_summary.csv"))
for _, r in ps.iterrows():
    expect(f"null mean {r.design}/{r.null}", f3(r["mean"]), MS + "\n" + SUPP)
c1w = ps[(ps.design == "C1") & (ps.null == "within_gene")].iloc[0]
expect("C1 within-gene p", f"p = {c1w.p_upper:.3f}")

ab = pd.read_csv(out("ablation_full.csv")).set_index("config")
for c in ("C1", "C8"):
    for m in ("lr", "rf", "gb"):
        expect(f"{c} {m}", f3(ab.loc[c, m]))
lr16 = pd.read_csv(out("ablation_lr16.csv")).set_index("config")
for c in lr16.index:
    expect(f"Table 1 {c}", f3(lr16.loc[c, "auc"]))

e = json.load(open(out("ensemble_audit.json")))
expect("ensemble reported", f3(e["auc_reported"]))
expect("ensemble nested", f3(e["auc_nested_blend"]))

wg = pd.read_csv(out("within_gene_discrimination.csv")).set_index("score")
for s in ("DCS", "RCS", "Conservation", "Normalised mobility", "AlphaMissense (variant)"):
    expect(f"within-gene {s}", f3(wg.loc[s, "auc_within"]))

tr = pd.read_csv(out("transfer_by_encoding.csv")).set_index("encoding")
for k in ("unrescaled", "conservation", "am_res", "am_res_pct"):
    expect(f"transfer {k}", f3(tr.loc[k, "auc_gene"]))
inc = pd.read_csv(out("incremental_over_am.csv")).set_index("model")
expect("AM alone gene-held-out", f3(inc.loc["AM alone", "auc"]))
expect("AM + published", f3(inc.loc["AM + position scores (published)", "auc"]))
gid = pd.read_csv(out("gene_identity_recovery.csv"))
r = gid[(gid.set == "labelled") & (gid.encoding == "published")].iloc[0]
expect("gene identity published", f"{r.bacc:.2f}")

# calibration, fold corrections and leave-two-genes-out (step 11)
BOTH = MS + "\n" + SUPP
cb = pd.read_csv(out("calibration_by_design.csv"))
for _, r in cb.iterrows():
    expect(f"calibration slope {r.model}/{r.split}", f"{r.cal_slope:+.2f}", BOTH)
fx = pd.read_csv(out("fold_corrections.csv"))
for _, r in fx.iterrows():
    expect(f"corrected AUC {r.model}/{r.correction}", f3(r.auc), BOTH)
lp = pd.read_csv(out("leave_two_genes_out.csv"))
for _, r in lp.iterrows():
    expect(f"leave-two-genes-out {r.model}", f3(r.auc_between_lpgo), BOTH)
cj = json.load(open(out("calibration.json")))
expect("fold correlation with training prevalence",
       f"{cj['fold_mean_pred_vs_train_prevalence_r']:.2f}", BOTH)

# influence, power and paired within-gene tests (step 12)
rj = json.load(open(out("robustness.json")))
expect("minimum detectable within-gene AUC", f"{rj['min_detectable_within_auc_80']:.2f}", BOTH)
pwr = pd.read_csv(out("within_gene_power.csv"))
for t in (0.65, 0.70):
    expect(f"power at {t}", f"{float(pwr[pwr.target_within_auc == t].power.iloc[0]):.2f}", BOTH)
pcm = pd.read_csv(out("within_gene_paired.csv"))
for _, r in pcm.iterrows():
    expect(f"paired within-gene {r.comparison[:28]}", f"{r.difference:+.3f}", BOTH)
gdi = pd.read_csv(out("gene_drop_influence.csv"))
expect("largest AUC after a single gene deletion", f3(gdi.pos_auc.max()), BOTH)
expect("smallest AUC after a single gene deletion", f3(gdi.pos_auc.min()), BOTH)

# a strictly increasing within-fold correction cannot change the within-gene component
offenders = [m for m in fx.model.unique()
             if fx.loc[fx.model == m, "auc_within"].round(4).nunique() != 1]
(ok if not offenders else bad).append(
    f"{'PASS' if not offenders else 'FAIL'}  fold corrections leave AUC_W unchanged"
    f"{'' if not offenders else ' (' + ', '.join(offenders) + ')'}")

# supplementary tables carry the same numbers
t1 = pd.read_csv(out("TableS1_discrimination.csv"))
for _, r in t1.iterrows():
    expect(f"Table S1 {r.score[:30]}", f"{r.auc:.3f}", SUPP)

# no unresolved placeholders or missing values
for name, text in (("manuscript", MS), ("supplement", SUPP)):
    for pat in [r"\bnan\b", r"\bNone\b", r"(?<![_^])\{[A-Za-z_]+\}", r"np\.float"]:
        hits = re.findall(pat, text)
        (ok if not hits else bad).append(
            f"{'PASS' if not hits else 'FAIL'}  {name}: no '{pat}' ({len(hits)} found)")

# citations: numbers appear in order of first use and every listed reference is cited
body = MS.split("# References")[0] if "# References" in MS else MS
nums = []
for m in re.findall(r"\[([0-9,–\-]+)\]", body):
    for part in m.split(","):
        part = part.replace("–", "-")
        if "-" in part:
            a, b = map(int, part.split("-"))
            nums.extend(range(a, b + 1))
        elif part.strip().isdigit():
            nums.append(int(part))
first = []
for n in nums:
    if n not in first:
        first.append(n)
listed = [int(x) for x in re.findall(r"^(\d+)\. ", MS, re.M)]
order_ok = first == sorted(first) and first == list(range(1, len(first) + 1))
(ok if order_ok else bad).append(f"{'PASS' if order_ok else 'FAIL'}  citations numbered in order of first use")
cover = set(first) == set(listed)
(ok if cover else bad).append(f"{'PASS' if cover else 'FAIL'}  every listed reference cited "
                              f"({len(listed)} listed, {len(set(first))} cited)")

print("=" * 72)
print("CONSISTENCY AUDIT")
print("=" * 72)
for line in bad:
    print(line)
print(f"\n{len(ok)} checks passed, {len(bad)} failed")
sys.exit(1 if bad else 0)
