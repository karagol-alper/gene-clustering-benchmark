#!/usr/bin/env python3
"""Step 10 - supplementary figures S1-S4."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from depbench import style as S
from depbench.common import FIGURES, out

S.apply()


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIGURES}/{name}.{ext}")
    plt.close(fig)


lab = pd.read_csv(out("cohort_full.csv"))
allv = pd.read_csv(out("clinvar_all_components.csv"))
maps = pd.read_csv(out("protein_scaling_maps.csv"))

# S1 ─ normalised mobility of ClinVar-observed residues per cohort protein
fig, ax = plt.subplots(figsize=(S.FULL_W, 2.9))
sub = allv[allv.Gene_Name.isin(lab.Gene_Name.unique())].drop_duplicates(["Protein_ID", "Residue_Pos"])
order = sub.groupby("Gene_Name").Norm_Dynamics.max().sort_values().index
rmsf_max = maps.set_index("Gene_Name").eval("rmsf_slope + rmsf_icpt")
rng = np.random.default_rng(0)
for i, gn in enumerate(order):
    v = sub.loc[sub.Gene_Name == gn, "Norm_Dynamics"].values
    ax.plot(i + rng.uniform(-0.18, 0.18, len(v)), v, "o", ms=2.2, color=S.BLUE, alpha=0.6, mew=0)
    lv = lab.loc[lab.Gene_Name == gn]
    ax.plot(np.full(len(lv), i + 0.32), lv.Norm_Dynamics, "_", color=S.ORANGE, ms=6, mew=1.2)
    ax.text(i, 1.04, f"{rmsf_max.get(gn, np.nan):.1f}", ha="center", fontsize=5.6, color=S.INK2)
ax.text(-0.9, 1.04, "max RMSF (Å):", ha="right", fontsize=5.6, color=S.INK2)
ax.set_xticks(range(len(order)), order, rotation=60, ha="right", fontsize=6.2)
ax.set_ylim(-0.03, 1.1)
ax.set_xlim(-1.0, len(order) - 0.4)
ax.set_ylabel("Normalised mobility (0-1 per protein)")
S.hgrid(ax)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=S.BLUE, ms=3, label="ClinVar-observed residues"),
                   Line2D([], [], marker="_", ls="", color=S.ORANGE, ms=6, mew=1.2, label="labelled variants")],
          loc="upper left", bbox_to_anchor=(0.27, 0.86), fontsize=6.0)
ax.set_title("Mobility is normalised by each protein's most mobile residue", loc="left", fontsize=8.5,
             fontweight="bold")
save(fig, "FigureS1")

# S2 ─ ablation for three model families with decomposition
ab = pd.read_csv(out("ablation_full.csv"))
fig, axes = plt.subplots(1, 3, figsize=(S.FULL_W, 2.5), sharey=True)
lbl = [f"{c}  {'V' if s == 'variant' else 'G'}|{'L' if a == 'full set' else 'F'}{'L' if b == 'full set' else 'F'}"
       for c, s, a, b in zip(ab.config, ab.split, ab.scaling, ab.selection)]
yy = np.arange(len(ab))[::-1]
for ax, (m, name) in zip(axes, [("lr", "Logistic regression (C = 1)"), ("rf", "Random forest"),
                               ("gb", "Gradient boosting")]):
    for yi, (_, r) in zip(yy, ab.iterrows()):
        c = S.ORANGE if r.split == "variant" else S.BLUE
        ax.plot([r[f"{m}_within"], r[m]], [yi, yi], color=S.GRID, lw=1.0, zorder=0)
        ax.plot(r[m], yi, "o", color=c, ms=4.2, mec="white", mew=0.5)
        ax.plot(r[f"{m}_within"], yi, "o", mfc="white", mec=c, ms=3.6, mew=0.9)
        ax.plot(r[f"{m}_between"], yi, "x", color=c, ms=3.4, mew=0.8)
    ax.axvline(0.5, color=S.BASE, lw=0.7)
    ax.set_xlim(0.45, 0.75)
    ax.set_title(name, loc="left", fontsize=7.5, fontweight="bold")
    ax.set_xlabel("AUC")
    S.vgrid(ax)
axes[0].set_yticks(yy, lbl, fontsize=6.0, family="DejaVu Sans Mono")
axes[2].legend(handles=[Line2D([], [], marker="o", ls="", color=S.INK2, label="pooled"),
                        Line2D([], [], marker="o", ls="", mfc="white", mec=S.INK2, label="within-gene"),
                        Line2D([], [], marker="x", ls="", color=S.INK2, label="between-gene")],
               loc="lower right", fontsize=5.8)
save(fig, "FigureS2")

# S3 ─ simulation: null case, logistic model and lookup
sim = pd.read_csv(out("simulation.csv"))
designs = ["12 genes x 2", "12 genes x 5", "12 genes x 10", "48 genes x 5", "cohort-matched (18 genes)"]
cols = [S.BLUE, S.ORANGE, S.AQUA, S.VIOLET, S.INK2]
fig, axes = plt.subplots(1, 3, figsize=(S.FULL_W, 2.45))
panels = [
    ("No individual signal (β = 0): local learner", 0.0, [("knn_vloo_mean", "-"), ("knn_gloo_mean", "--")], "AUC"),
    ("Logistic model on one shifted feature", 0.75, [("feat_infl_vloo_new_mean", "-"), ("feat_bias_gloo_new_mean", "--")], "Bias relative to new-gene AUC"),
    ("Gene-prevalence lookup (β = 0)", 0.0, [("prev_vloo_mean", "-"), ("prev_gloo_mean", "--")], "AUC"),
]
for ax, (title, beta, series, ylab) in zip(axes, panels):
    sel = sim[(sim.beta == beta) & (sim.gamma == 1.0)]
    for d_, c in zip(designs, cols):
        dd = sel[sel.design == d_].sort_values("gene_sd")
        for col, ls in series:
            ax.plot(dd.gene_sd, dd[col], ls, color=c, lw=1.2, marker="o", ms=2.6,
                    label=d_.replace(" x ", " × ") if ls == "-" else None)
    ax.axhline(0.5 if ylab == "AUC" else 0.0, color=S.BASE, lw=0.8)
    ax.set_title(title, loc="left", fontsize=7.2, fontweight="bold")
    ax.set_xlabel("Between-gene SD of log-odds")
    ax.set_ylabel(ylab)
    S.hgrid(ax)
handles, labels = axes[0].get_legend_handles_labels()
handles += [Line2D([], [], color=S.INK2, ls="-", lw=1.2), Line2D([], [], color=S.INK2, ls="--", lw=1.2)]
labels += ["variant-level CV", "gene-level CV"]
fig.legend(handles, labels, loc="lower center", ncol=7, fontsize=5.8, bbox_to_anchor=(0.5, -0.06),
           columnspacing=1.0, handlelength=1.8)
fig.tight_layout(rect=(0, 0.06, 1, 1))
save(fig, "FigureS3")

# S4 ─ cluster versus naive bootstrap width
t1 = pd.read_csv(out("TableS1_discrimination.csv"))
fig, ax = plt.subplots(figsize=(S.FULL_W * 0.62, 3.0))
t1 = t1.iloc[::-1].reset_index(drop=True)
ax.barh(range(len(t1)), t1.width_ratio, color=[S.BLUE if "†" in s else S.ORANGE for s in t1.score],
        height=0.66)
ax.axvline(1, color=S.INK2, lw=0.8)
ax.set_yticks(range(len(t1)), t1.score, fontsize=5.9)
ax.set_xlabel("Width of gene-cluster / naive bootstrap 95% interval")
S.vgrid(ax)
ax.legend(handles=[Line2D([], [], color=S.ORANGE, lw=5, label="raw score or reported prediction"),
                   Line2D([], [], color=S.BLUE, lw=5, label="gene-held-out model")],
          fontsize=5.8, loc="lower right")
save(fig, "FigureS4")
print("supplementary figures written")
