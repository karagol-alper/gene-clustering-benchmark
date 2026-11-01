#!/usr/bin/env python3
"""Step 13 - Figure 6: cross-fold calibration and what the cohort can resolve."""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from depbench import style as S
from depbench.common import FIGURES, out

S.apply()
fcal = pd.read_csv(out("fold_calibration.csv"))
fcal = fcal[fcal.model == "Position-level scores"]
fx = pd.read_csv(out("fold_corrections.csv"))
lp = pd.read_csv(out("leave_two_genes_out.csv")).set_index("model")
pwr = pd.read_csv(out("within_gene_power.csv"))
pcm = pd.read_csv(out("within_gene_paired.csv"))
rj = json.load(open(out("robustness.json")))
wg = pd.read_csv(out("within_gene_discrimination.csv")).set_index("score")

fig = plt.figure(figsize=(S.FULL_W, 2.75))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.08, 1.0], wspace=0.46)

# ── (a) each fold's mean prediction against its training prevalence ─────
ax = fig.add_subplot(gs[0])
S.hgrid(ax)
pure_p = fcal.gene_prevalence == 1
pure_b = fcal.gene_prevalence == 0
ax.scatter(fcal.train_prevalence[pure_b], fcal.mean_prediction[pure_b], s=26, color=S.BEN,
           edgecolor="white", lw=0.7, zorder=3, label="held-out gene all B/LB")
ax.scatter(fcal.train_prevalence[pure_p], fcal.mean_prediction[pure_p], s=26, color=S.PATH,
           edgecolor="white", lw=0.7, zorder=3, label="held-out gene all P/LP")
mix = ~(pure_p | pure_b)
ax.scatter(fcal.train_prevalence[mix], fcal.mean_prediction[mix], s=26, color=S.MUTED,
           edgecolor="white", lw=0.7, zorder=3, label="both classes")
b, a = np.polyfit(fcal.train_prevalence, fcal.mean_prediction, 1)
xx = np.linspace(fcal.train_prevalence.min(), fcal.train_prevalence.max(), 20)
ax.plot(xx, a + b * xx, color=S.INK2, lw=1.0, ls="--", zorder=2)
r = np.corrcoef(fcal.train_prevalence, fcal.mean_prediction)[0, 1]
ax.annotate(f"r = {r:.2f}", xy=(0.97, 0.06), xycoords="axes fraction", ha="right",
            fontsize=7, color=S.INK2)
hm = fcal[fcal.gene == "HMGCS2"]
if len(hm):
    ax.annotate("HMGCS2", xy=(hm.train_prevalence.iloc[0], hm.mean_prediction.iloc[0]),
                xytext=(6, -2), textcoords="offset points", fontsize=6.4, color=S.PATH)
ax.set_xlabel("pathogenic fraction of the training genes")
ax.set_ylabel("mean prediction for the held-out gene")
ax.legend(loc="upper left", fontsize=6.2, handletextpad=0.3, borderpad=0.2)
S.panel(ax, "a", "Folds track their prior")

# ── (b) corrections that do not compare across folds ────────────────────
ax = fig.add_subplot(gs[1])
S.vgrid(ax)
models = ["Position-level scores", "Full 19-feature model", "AlphaMissense only"]
short = ["Position-level\nscores", "Full 19-feature\nmodel", "AlphaMissense\nonly"]
labels = ["pooled LOGO", "prior shift", "leave-two-genes-out"]
cols = [S.MUTED, S.BLUE, S.AQUA]
h = 0.24
for i, m in enumerate(models):
    base = fx[(fx.model == m) & (fx.correction == "none")].iloc[0]
    pri = fx[(fx.model == m) & (fx.correction == "prior shift")].iloc[0]
    vals = [base.auc_between, pri.auc_between, lp.loc[m, "auc_between_lpgo"]]
    for j, v in enumerate(vals):
        ax.barh(i + (1 - j) * h, v, height=h * 0.9, color=cols[j], edgecolor="white", lw=0.6,
                zorder=3)
ax.axvline(0.5, color=S.INK2, lw=0.9, ls=":", zorder=2)
ax.set_yticks(range(len(models)), short, fontsize=6.6)
ax.set_xlabel("between-gene component, AUC$_B$")
ax.set_xlim(0, 1)
ax.invert_yaxis()
ax.legend(handles=[Line2D([], [], color=c, lw=4, label=l) for c, l in zip(cols, labels)],
          loc="upper right", fontsize=6.2, handlelength=1.4, handletextpad=0.4,
          borderpad=0.25, labelspacing=0.35)
S.panel(ax, "b", "Cross-fold comparison removed")

# ── (c) power of the within-gene test in this cohort ────────────────────
ax = fig.add_subplot(gs[2])
S.hgrid(ax)
ax.plot(pwr.target_within_auc, pwr.power, color=S.BLUE, marker="o", ms=3.4, zorder=3)
ax.axhline(0.8, color=S.INK2, lw=0.9, ls=":", zorder=2)
mde = rj["min_detectable_within_auc_80"]
ax.plot([mde, mde], [0, 0.8], color=S.INK2, lw=0.9, ls=":", zorder=2)
ax.annotate(f"80% power at\nAUC$_W$ = {mde:.2f}", xy=(0.97, 0.22), xycoords="axes fraction",
            ha="right", va="center", fontsize=6.4, color=S.INK2)
marks = [("conservation", "Conservation", S.VIOLET, 1.03, "right", -2),
         ("DCS", "DCS", S.ORANGE, 1.11, "left", 2),
         ("AlphaMissense", "AlphaMissense (variant)", S.AQUA, 1.03, "right", -2)]
for name, key, col, ytop, ha, dx in marks:
    v = float(wg.loc[key, "auc_within"])
    ax.plot([v, v], [0, ytop], color=col, lw=1.1, zorder=2)
    ax.annotate(name, xy=(v, ytop), xytext=(dx, 1), textcoords="offset points", ha=ha,
                fontsize=6.1, color=col)
ax.set_xlabel("true within-gene AUC")
ax.set_ylabel("power of the stratified permutation test")
ax.set_ylim(0, 1.22)
ax.set_xlim(0.49, 0.945)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
S.panel(ax, "c", "Detectable effects")

for ext in ("png", "pdf"):
    fig.savefig(f"{FIGURES}/Figure6.{ext}")
fig.savefig(f"{FIGURES}/Figure6.tiff", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print("Figure 6 written")
