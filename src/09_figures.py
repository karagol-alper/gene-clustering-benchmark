#!/usr/bin/env python3
"""Step 9 - main-text figures (600 dpi PNG, vector PDF and 300 dpi TIFF)."""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.metrics import roc_curve

from depbench import style as S
from depbench.common import FIGURES, out

S.apply()
lab = pd.read_csv(out("cohort_full.csv"))
y, g = lab.label.values, lab.Gene_Name.values
diag = json.load(open(out("diagnostics.json")))
ens = json.load(open(out("ensemble_audit.json")))
V_PAIRS = diag["pairs"]["n_within"]


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIGURES}/{name}.{ext}")
    fig.savefig(f"{FIGURES}/{name}.tiff", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


# ═════════════════════════════════════════════════════════════ Figure 1
def figure1():
    fig = plt.figure(figsize=(S.FULL_W, 2.65))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.82, 1.18], wspace=0.55)
    # (a) composition
    ax = fig.add_subplot(gs[0])
    comp = lab.groupby("Gene_Name").label.agg(["size", "sum"]).sort_values(["size", "sum"])
    yy = np.arange(len(comp))
    ax.barh(yy, comp["sum"], color=S.PATH, height=0.72, label="P/LP", edgecolor="white", lw=0.6)
    ax.barh(yy, comp["size"] - comp["sum"], left=comp["sum"], color=S.BEN, height=0.72,
            label="B/LB", edgecolor="white", lw=0.6)
    ax.set_yticks(yy, comp.index, fontsize=6.3)
    for i, (gn, r) in enumerate(comp.iterrows()):
        mixed = 0 < r["sum"] < r["size"]
        ax.text(r["size"] + 0.6, i, "mixed" if mixed else "", va="center", fontsize=5.8, color=S.INK2)
    ax.set_xlabel("Variants")
    ax.set_xlim(0, 38)
    S.vgrid(ax)
    ax.legend(loc="lower right", handlelength=1.0)
    S.panel(ax, "a", "Variants per gene", x=-0.3)

    # (b) pair composition
    ax = fig.add_subplot(gs[1])
    pr = diag["pairs"]
    per = pr["per_gene"]
    segs = [("Between genes", pr["n_between"], "#c9c8c1")] + \
           [(f"Within {k}", v, c) for (k, v), c in zip(sorted(per.items(), key=lambda t: -t[1]),
                                                     [S.VIOLET, S.AQUA, S.YELLOW])]
    bottom = 0
    for lbl, n, c in segs:
        ax.bar(0, n, bottom=bottom, color=c, width=0.55, edgecolor="white", lw=0.8)
        bottom += n
    ax.set_xlim(-0.5, 1.9)
    ax.set_ylim(0, pr["n_pairs"] * 1.02)
    ax.set_xticks([])
    ax.set_ylabel("Pathogenic-benign pairs")
    ax.spines["bottom"].set_visible(False)
    ax.text(0.34, pr["n_between"] / 2, f"between genes\n{pr['n_between']:,} "
            f"({100 * (1 - pr['w_within']):.0f}%)", va="center", fontsize=6.5, color=S.INK)
    ax.text(0.34, pr["n_between"] + pr["n_within"] / 2,
            f"within genes\n{pr['n_within']} ({100 * pr['w_within']:.0f}%)\n" +
            ", ".join(f"{k} {v}" for k, v in sorted(per.items(), key=lambda t: -t[1])),
            va="center", fontsize=6.0, color=S.INK)
    S.panel(ax, "b", "What the AUC compares", x=-0.55)

    # (c) within vs between
    ax = fig.add_subplot(gs[2])
    wg = pd.read_csv(out("within_gene_discrimination.csv")).set_index("score")
    pts = [("AlphaMissense (variant)", "AM variant", S.BLUE, (5, -3), True),
           ("AlphaMissense (residue mean)", "AM residue", S.BLUE, (5, -10), True),
           ("DCS", "DCS", S.ORANGE, (5, 2), True), ("RCS", "RCS", S.ORANGE, (5, -9), True),
           ("Conservation", "conservation", S.ORANGE, (-7, 0), True),
           ("Normalised mobility", "mobility", S.ORANGE, (5, -3), True),
           ("|delta hydrophobicity|", "|Δhydro.|", S.MUTED, (0, -12), False)]
    ax.axhline(0.5, color=S.BASE, lw=0.7, zorder=0)
    ax.axvline(0.5, color=S.BASE, lw=0.7, zorder=0)
    ax.plot([0, 1], [0, 1], color=S.GRID, lw=0.8, zorder=0)
    for key, lbl, c, off, bars in pts:
        r = wg.loc[key]
        lo, hi = map(float, r.within_ci.split("-"))
        if bars:
            ax.errorbar(r.auc_between, r.auc_within, yerr=[[r.auc_within - lo], [hi - r.auc_within]],
                        fmt="o", color=c, ms=4.2, elinewidth=0.8, capsize=0, mec="white", mew=0.6)
        else:
            ax.plot(r.auc_between, r.auc_within, "o", color=c, ms=4.2, mec="white", mew=0.6)
        ax.annotate(lbl, (r.auc_between, r.auc_within), xytext=off, textcoords="offset points",
                    fontsize=6.2, color=S.INK2,
                    ha="right" if off[0] < 0 else ("center" if off[0] == 0 else "left"))
    prev = diag["prevalence"]
    for key, lbl, off in [("in_sample", "lookup, in sample", (-6, -11)),
                          ("variant_loo", "lookup, variant excluded", (-7, 2))]:
        dd = prev[key]
        ax.plot(dd["auc_between"], dd["auc_within"], "D", color=S.AQUA, ms=4.4, mec="white", mew=0.6)
        ax.annotate(lbl, (dd["auc_between"], dd["auc_within"]), xytext=off, textcoords="offset points",
                    fontsize=6.2, color=S.INK2, ha="right")
    e = ens["ensemble_decomposition"]
    ax.plot(e["auc_between"], e["auc_within"], "s", color=S.AQUA, ms=4.4, mec="white", mew=0.6)
    ax.annotate("tuned ensemble", (e["auc_between"], e["auc_within"]), xytext=(4, -26),
                textcoords="offset points", fontsize=6.2, color=S.INK2, ha="left",
                arrowprops=dict(arrowstyle="-", color=S.MUTED, lw=0.5, shrinkA=0, shrinkB=3))
    ax.set_xlim(0.3, 1.0)
    ax.set_ylim(-0.04, 1.02)
    ax.set_xlabel("Between-gene AUC")
    ax.set_ylabel("Within-gene AUC")
    S.panel(ax, "c", "Two components of one AUC", x=-0.2)
    save(fig, "Figure1")


# ═════════════════════════════════════════════════════════════ Figure 2
def figure2():
    fig = plt.figure(figsize=(S.FULL_W, 3.05))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.2, 1], wspace=0.55, left=0.06, right=0.99,
                          top=0.88, bottom=0.27)
    from depbench.metrics import fast_auc
    # (a) ROC
    ax = fig.add_subplot(gs[0])
    gp = pd.read_csv(out("gene_prevalence_scores.csv"))
    ep = pd.read_csv(out("ensemble_predictions.csv"))
    curves = [(gp.in_sample, gp.label, "Lookup, in sample", S.AQUA, "-"),
              (gp.variant_loo, gp.label, "Lookup, variant excluded", S.AQUA, "--"),
              (ep.ensemble_reported, ep.label, "Ensemble, as reported", S.ORANGE, "-"),
              (ep.ensemble_nested_blend, ep.label, "Ensemble, nested blend", S.ORANGE, ":")]
    for sc, yy, lbl, c, ls in curves:
        f, t, _ = roc_curve(yy, sc)
        ax.plot(f, t, color=c, ls=ls, lw=1.4, label=f"{lbl} ({fast_auc(yy, sc):.3f})")
    ax.plot([0, 1], [0, 1], color=S.BASE, lw=0.7)
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.2), fontsize=5.8, handlelength=2.0, ncol=1)
    S.panel(ax, "a", "Lookup versus ensemble", x=-0.18)

    # (b) 16 configurations
    ax = fig.add_subplot(gs[1])
    lr = pd.read_csv(out("ablation_lr16.csv"))
    lr["lbl"] = [f"{'V' if a == 'variant' else 'G'} | {'L' if b == 'full set' else 'F'}"
                 f"{'L' if c == 'full set' else 'F'}{'L' if d == 'full set' else 'F'}"
                 for a, b, c, d in zip(lr.split, lr.scaling, lr.selection, lr.tuning)]
    lr = lr.iloc[::-1].reset_index(drop=True)
    yy = np.arange(len(lr))
    for i, r in lr.iterrows():
        c = S.ORANGE if r.split == "variant" else S.BLUE
        ax.plot([r.auc_within, r.auc], [i, i], color=S.GRID, lw=1.0, zorder=0)
        ax.plot(r.auc, i, "o", color=c, ms=4.2, mec="white", mew=0.5)
        ax.plot(r.auc_within, i, "o", mfc="white", mec=c, ms=3.4, mew=0.9)
    ax.set_yticks(yy, lr.lbl, fontsize=5.8, family="DejaVu Sans Mono")
    ax.set_xlim(0.5, 0.72)
    ax.set_ylim(-0.7, len(lr) - 0.3)
    ax.set_xlabel("AUC (logistic regression)")
    S.vgrid(ax)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=S.ORANGE, label="variant-level split"),
                       Line2D([], [], marker="o", ls="", color=S.BLUE, label="gene-level split"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec=S.INK2, label="within-gene component")],
              loc="upper left", bbox_to_anchor=(-0.05, -0.2), fontsize=5.8, ncol=2)
    ax.text(1.0, -0.42, "split | scaling, selection, tuning\nL = full labelled set, F = in fold",
            transform=ax.transAxes, fontsize=5.6, color=S.INK2, va="top", ha="right")
    S.panel(ax, "b", "Factorial ablation", x=-0.25)

    # (c) permutation nulls, variant-level C1
    ax = fig.add_subplot(gs[2])
    dist = pd.read_csv(out("permutation_distributions.csv"))
    summ = pd.read_csv(out("permutation_summary.csv"))
    obs = float(summ.loc[summ.design == "C1", "observed"].iloc[0])
    bins = np.linspace(0.25, 0.9, 40)
    for col, lbl, c in [("C1_global", "labels permuted globally", S.MUTED),
                        ("C1_within_gene", "labels permuted within gene", S.VIOLET)]:
        v = dist[col].dropna()
        ax.hist(v, bins=bins, color=c, alpha=0.55, label=lbl, edgecolor="white", lw=0.3)
    ax.axvline(obs, color=S.ORANGE, lw=1.4)
    ax.axvline(0.5, color=S.BASE, lw=0.7, ls=":")
    pg = summ[(summ.design == "C1") & (summ.null == "global")].iloc[0]
    pw = summ[(summ.design == "C1") & (summ.null == "within_gene")].iloc[0]
    top = ax.get_ylim()[1]
    ax.set_ylim(0, top * 1.05)
    ax.set_xlim(0.25, 0.97)
    ax.text(obs + 0.012, top * 0.99, f"observed {obs:.3f}", fontsize=6.0, color=S.INK, va="top")
    ax.text(obs + 0.012, top * 0.86, f"p (global)\n= {pg.p_upper:.3f}\n\np (within gene)\n= {pw.p_upper:.3f}",
            fontsize=5.8, color=S.INK2, va="top")
    ax.set_xlabel("Variant-level AUC, configuration C1")
    ax.set_ylabel("Replicates")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.2), fontsize=5.8)
    S.panel(ax, "c", "Permutation nulls", x=-0.2)
    save(fig, "Figure2")


# ═════════════════════════════════════════════════════════════ Figure 3
def figure3():
    fig = plt.figure(figsize=(S.FULL_W, 3.15))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.9, 1.15, 1.15], wspace=0.95, left=0.06, right=0.95,
                          top=0.89, bottom=0.27)
    allv = pd.read_csv(out("clinvar_all_components.csv"))
    ax = fig.add_subplot(gs[0])
    sub = allv[allv.Gene_Name.isin(lab.Gene_Name.unique())].drop_duplicates(["Protein_ID", "Residue_Pos"])
    mult = sub.groupby("Gene_Name").dcs_scale_factor.first().sort_values()
    hi = [mult.index[0], "PGK1", mult.index[-1]]
    for gn, d in sub.groupby("Gene_Name"):
        if gn in hi:
            continue
        ax.plot(d.Raw_DCS, d.DCS_Score, "o", ms=1.6, color="#c9c8c1", mew=0)
    for gn, c in zip(hi, [S.BLUE, S.AQUA, S.ORANGE]):
        d = sub[sub.Gene_Name == gn]
        ax.plot(d.Raw_DCS, d.DCS_Score, "o", ms=2.3, color=c, mew=0, label=f"{gn} (×{mult[gn]:.1f})")
    ax.set_xlabel("Conservation × normalised mobility")
    ax.set_ylabel("DCS as published (0–10)")
    ax.set_xlim(0, 0.45)
    ax.legend(loc="upper left", bbox_to_anchor=(-0.08, -0.2), fontsize=5.8, handletextpad=0.2,
              title="protein (multiplier); grey, other proteins", title_fontsize=5.6)
    S.panel(ax, "a", "Per-protein rescaling", x=-0.3)

    ax = fig.add_subplot(gs[1])
    gid = pd.read_csv(out("gene_identity_recovery.csv"))
    order = ["published", "unrescaled", "components", "conservation", "mobility", "rmsf",
             "pct_protein", "am_res", "am_res_minmax", "am_res_xdcs", "am_res_pct", "substitution"]
    names = {"published": "published scores", "unrescaled": "before rescaling",
             "components": "conservation + mobility", "conservation": "conservation",
             "mobility": "normalised mobility", "rmsf": "mean RMSF",
             "pct_protein": "within-protein percentile", "am_res": "AM residue mean",
             "am_res_minmax": "AM residue, min–max", "am_res_xdcs": "AM residue × DCS mult.",
             "am_res_pct": "AM residue, percentile", "substitution": "substitution descriptors"}
    L = gid[gid.set == "labelled"].set_index("encoding")
    A = gid[gid.set == "all_clinvar"].set_index("encoding")
    yy = np.arange(len(order))[::-1]
    for yi, e in zip(yy, order):
        c = S.ORANGE if e in order[:7] else (S.BLUE if e.startswith("am") else S.MUTED)
        r = L.loc[e]
        ax.plot([r.null_mean - 1.96 * r.null_sd, r.null_mean + 1.96 * r.null_sd], [yi, yi],
                color=S.GRID, lw=3.2, solid_capstyle="butt", zorder=0)
        ax.plot(r.bacc, yi, "o", color=c, ms=4.2, mec="white", mew=0.5)
        if e in A.index:
            ax.plot(A.loc[e].bacc, yi, "o", mfc="white", mec=c, ms=3.4, mew=0.9)
    ax.set_yticks(yy, [names[e] for e in order], fontsize=6.0)
    ax.axvline(1 / 8, color=S.BASE, lw=0.7, ls=":")
    ax.set_xlabel("Balanced accuracy for gene identity")
    ax.set_xlim(-0.01, 0.3)
    S.vgrid(ax)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=S.INK2, label="labelled cohort (8 genes)"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec=S.INK2,
                              label=f"all ClinVar residues ({int(A.k.iloc[0])} proteins)"),
                       Line2D([], [], color=S.GRID, lw=3.2, label="permutation null, 95% (labelled)")],
              loc="upper left", bbox_to_anchor=(-0.1, -0.2), fontsize=5.6)
    S.panel(ax, "b", "Protein identity", x=-0.6)

    ax = fig.add_subplot(gs[2])
    wg = pd.read_csv(out("within_gene_discrimination.csv")).set_index("score")
    keys = [("DCS", "DCS", S.ORANGE), ("RCS", "RCS", S.ORANGE), ("Conservation", "conservation", S.ORANGE),
            ("Normalised mobility", "normalised mobility", S.ORANGE), ("Mean RMSF", "mean RMSF", S.ORANGE),
            ("|delta hydrophobicity|", "|Δhydrophobicity|", S.MUTED),
            ("AlphaMissense (residue mean)", "AM residue mean", S.BLUE),
            ("AlphaMissense (variant)", "AM variant", S.BLUE)]
    yy = np.arange(len(keys))[::-1]
    for yi, (k, lbl, c) in zip(yy, keys):
        r = wg.loc[k]
        lo, hi_ = map(float, r.within_ci.split("-"))
        ax.plot([lo, hi_], [yi, yi], color=c, lw=1.0, alpha=0.6)
        ax.plot(r.auc_within, yi, "o", color=c, ms=4.4, mec="white", mew=0.5)
        for gn, mk in [("PGK1", "|"), ("GHRHR", "x"), ("KAT8", "+")]:
            v = r.get(f"auc_{gn}")
            if pd.notna(v):
                ax.plot(v, yi + 0.3, mk, color=S.INK2, ms=3.0, mew=0.7)
        pv = r.within_p_two_sided
        ax.text(1.04, yi, f"p = {pv:.3f}" if pv >= 0.001 else "p < 0.001", fontsize=5.8, va="center",
                color=S.INK2, transform=ax.get_yaxis_transform())
    ax.set_yticks(yy, [k[1] for k in keys], fontsize=6.0)
    ax.axvline(0.5, color=S.BASE, lw=0.7)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.6, len(keys) - 0.3)
    ax.set_xlabel(f"Within-gene AUC ({V_PAIRS} pairs)")
    ax.legend(handles=[Line2D([], [], marker=m, ls="", color=S.INK2, label=n, ms=3.2)
                       for m, n in [("|", "PGK1 (288 pairs)"), ("x", "GHRHR (6)"), ("+", "KAT8 (4)")]],
              loc="upper left", bbox_to_anchor=(-0.1, -0.2), fontsize=5.6, ncol=3, columnspacing=0.8)
    S.panel(ax, "c", "Signal within genes", x=-0.55)
    save(fig, "Figure3")


# ═════════════════════════════════════════════════════════════ Figure 4
def figure4():
    fig = plt.figure(figsize=(S.FULL_W, 3.1))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 0.75, 1.0], wspace=1.05, left=0.165, right=0.93,
                          top=0.89, bottom=0.24)
    ax = fig.add_subplot(gs[0])
    tr = pd.read_csv(out("transfer_by_encoding.csv")).set_index("encoding")
    order = [("published", "published scores"), ("unrescaled", "before rescaling"),
             ("components", "conservation + mobility"), ("conservation", "conservation"),
             ("mobility", "normalised mobility"), ("pct_protein", "within-protein percentile"),
             ("pct_cohort", "within-gene percentile †"), ("z_cohort", "within-gene z-score †"),
             ("am_res", "AM residue mean"), ("am_res_minmax", "AM residue, min–max"),
             ("am_res_xdcs", "AM residue × DCS mult."), ("am_res_pct", "AM residue, percentile")]
    yy = np.arange(len(order))[::-1]
    for yi, (e, lbl) in zip(yy, order):
        r = tr.loc[e]
        c = S.BLUE if e.startswith("am") else S.ORANGE
        lo_, hi_ = min(r.within_gene, r.between_gene, r.auc_gene), max(r.within_gene, r.between_gene, r.auc_gene)
        ax.plot([lo_, hi_], [yi, yi], color=S.GRID, lw=1.0, zorder=0)
        ax.plot(r.auc_gene, yi, "o", color=c, ms=4.4, mec="white", mew=0.5, zorder=3)
        ax.plot(r.within_gene, yi + 0.18, "^", mfc="white", mec=c, ms=3.4, mew=0.8)
        ax.plot(r.between_gene, yi - 0.18, "v", mfc="white", mec=c, ms=3.4, mew=0.8)
    ax.set_yticks(yy, [o[1] for o in order], fontsize=6.0)
    ax.axvline(0.5, color=S.BASE, lw=0.7)
    ax.set_xlim(0.05, 0.95)
    ax.set_xlabel("Gene-held-out AUC")
    S.vgrid(ax)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=S.INK2, label="pooled"),
                       Line2D([], [], marker="^", ls="", mfc="white", mec=S.INK2, label="within-gene pairs"),
                       Line2D([], [], marker="v", ls="", mfc="white", mec=S.INK2, label="between-gene pairs")],
              loc="upper left", bbox_to_anchor=(-0.05, -0.17), fontsize=5.6, ncol=3, columnspacing=0.8)
    S.panel(ax, "a", "Transfer to held-out genes", x=-0.62)

    ax = fig.add_subplot(gs[1])
    dist = pd.read_csv(out("permutation_distributions.csv"))
    summ = pd.read_csv(out("permutation_summary.csv"))
    obs = float(summ.loc[summ.design == "POS", "observed"].iloc[0])
    data = [dist[c].dropna().values for c in ["POS_global", "POS_donor", "POS_within_gene"]]
    parts = ax.violinplot(data, positions=[0, 1, 2], widths=0.78, showextrema=False)
    for b, c in zip(parts["bodies"], [S.MUTED, S.AQUA, S.VIOLET]):
        b.set_facecolor(c); b.set_alpha(0.5); b.set_edgecolor("none")
    for i, d in enumerate(data):
        ax.plot([i - 0.2, i + 0.2], [np.mean(d)] * 2, color=S.INK, lw=1.0)
    ax.axhline(obs, color=S.ORANGE, lw=1.3)
    floor = diag["prevalence"]["gene_loo"]["auc"]
    ax.axhline(floor, color=S.AQUA, lw=1.0, ls="--")
    ax.axhline(0.5, color=S.BASE, lw=0.7, ls=":")
    ax.text(-0.55, obs - 0.012, f"observed {obs:.3f}", fontsize=5.8, va="top", color=S.INK)
    ax.text(-0.55, floor + 0.012, f"prevalence only {floor:.3f}", fontsize=5.8, va="bottom",
            color=S.INK2)
    ax.set_xticks([0, 1, 2], ["global", "donor", "within\ngene"], fontsize=6.0)
    ax.set_xlabel("Null")
    ax.set_xlim(-0.6, 2.5)
    ax.set_ylabel("Gene-held-out AUC, position scores")
    ax.set_ylim(0, 0.85)
    S.panel(ax, "b", "Null calibration", x=-0.38)

    ax = fig.add_subplot(gs[2])
    inc = pd.read_csv(out("incremental_over_am.csv"))
    inc = inc[inc.model != "AM alone"].copy()
    lbls = {"AM + position scores (published)": "+ published scores",
            "AM + position scores (unrescaled)": "+ before rescaling",
            "AM + conservation": "+ conservation",
            "AM + position scores (within-protein percentile)": "+ protein percentile",
            "AM + position scores (within-gene percentile, labelled)": "+ gene percentile †",
            "AM + position scores (within-gene z, labelled)": "+ gene z-score †",
            "AM + position scores (within-gene centred, labelled)": "+ gene centred †",
            "AM + substitution descriptors": "+ substitutions",
            "AM + full feature set": "+ full feature set"}
    inc = inc[inc.model.isin(lbls)]
    inc["lbl"] = inc.model.map(lbls)
    inc = inc.iloc[::-1].reset_index(drop=True)
    for i, r in inc.iterrows():
        c = S.ORANGE if ("position" in r.model or "conservation" in r.model) else S.MUTED
        ax.plot([r.delta_cluster_ci_lo, r.delta_cluster_ci_hi], [i, i], color=c, lw=1.0, alpha=0.7)
        ax.plot(r.delta, i, "o", color=c, ms=4.2, mec="white", mew=0.5)
        ax.text(1.04, i, f"p = {r.clustered_p:.2f}", fontsize=5.8, va="center", color=S.INK2,
                transform=ax.get_yaxis_transform())
    ax.set_yticks(range(len(inc)), inc.lbl, fontsize=6.0)
    ax.axvline(0, color=S.BASE, lw=0.8)
    base = pd.read_csv(out("incremental_over_am.csv")).set_index("model").loc["AM alone", "auc"]
    ax.set_xlabel(f"Change in gene-held-out AUC\n(AlphaMissense alone {base:.3f})")
    S.vgrid(ax)
    S.panel(ax, "c", "Increment over AM", x=-0.62)
    save(fig, "Figure4")


# ═════════════════════════════════════════════════════════════ Figure 5
def figure5():
    sim = pd.read_csv(out("simulation.csv"))
    req = pd.read_csv(out("design_requirements.csv"))
    fig = plt.figure(figsize=(S.FULL_W, 2.75))
    gs = fig.add_gridspec(1, 3, wspace=0.45, left=0.07, right=0.99, top=0.88, bottom=0.17)
    designs = ["12 genes x 5", "12 genes x 10", "48 genes x 5", "cohort-matched (18 genes)"]
    cols = [S.BLUE, S.ORANGE, S.AQUA, S.INK2]
    sel = sim[(sim.beta == 0.75) & (sim.gamma == 1.0)]
    for k, (metric, title, ylab) in enumerate([
            ("knn_infl_vloo_new_mean", "Variant-level CV", "Variant-level AUC − new-gene AUC"),
            ("knn_bias_gloo_new_mean", "Gene-level CV", "Gene-held-out AUC − new-gene AUC")]):
        ax = fig.add_subplot(gs[k])
        for d_, c in zip(designs, cols):
            dd = sel[sel.design == d_].sort_values("gene_sd")
            ax.plot(dd.gene_sd, dd[metric], "-o", color=c, ms=3.0, lw=1.3,
                    label=d_.replace(" x ", " × ").replace("cohort-matched (18 genes)", "observed gene sizes (18)"))
        if k == 1:
            w = sel[sel.design.isin(designs)].groupby("gene_sd").knn_bias_gloo_w_mean.mean()
            ax.plot(w.index, w.values, ":", color=S.VIOLET, lw=1.4, label="within-gene component\n(mean over designs)")
            ax.legend(loc="lower left", fontsize=5.6, handlelength=1.6)
        ax.axhline(0, color=S.BASE, lw=0.8)
        ax.set_xlabel("Between-gene SD of log-odds")
        ax.set_ylabel(ylab)
        S.hgrid(ax)
        if k == 0:
            ax.legend(loc="upper left", fontsize=5.6)
        S.panel(ax, "ab"[k], title, x=-0.2)
    ax = fig.add_subplot(gs[2])
    for est, c, nm in [("within", S.VIOLET, "within-gene AUC"), ("pooled", S.MUTED, "pooled AUC")]:
        for h, mk, ls in [(0.05, "o", "-"), (0.03, "s", "--")]:
            a = req[(req.estimator == est) & (req.target_within_auc == 0.65) & (req.half_width == h)].sort_values("gene_sd")
            ax.plot(a.gene_sd, a.genes_required, ls, marker=mk, color=c, ms=3.4, lw=1.3, label=f"{nm}, ±{h}")
    ax.axhline(18, color=S.ORANGE, lw=1.0)
    ax.text(1.0, 16.0, "this cohort: 18 genes, 3 with both classes", fontsize=5.8, color=S.INK2, va="top")
    ax.set_yscale("log")
    ax.set_ylim(10, 2500)
    ax.set_xticks([1, 2, 4])
    ax.set_xlabel("Between-gene SD of log-odds")
    ax.set_ylabel("Genes required (5 variants each)")
    ax.legend(loc="lower right", fontsize=5.5, bbox_to_anchor=(1.0, 0.12))
    S.hgrid(ax)
    S.panel(ax, "c", "Genes needed", x=-0.2)
    save(fig, "Figure5")


if __name__ == "__main__":
    import sys
    which = sys.argv[1:] or ["1", "2", "3", "4", "5"]
    for w in which:
        globals()[f"figure{w}"]()
        print(f"Figure {w} written")
