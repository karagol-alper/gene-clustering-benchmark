#!/usr/bin/env python3
"""
Step 0 - build the labelled cohort and reconstruct the components of the
position-level scores.

The published position-level scores are derived quantities:
    Raw_DCS        = conservation x Norm_Dynamics
    Raw_RCS        = conservation x (1 - Norm_Dynamics)
    DCS_Score      = per-protein min-max rescaling of Raw_DCS to 0-10
    RCS_Score      = per-protein min-max rescaling of Raw_RCS to 0-10
    Normalized_DCS = log(1 + DCS_Score) / log(11)
where Norm_Dynamics is itself per-protein-normalised mean backbone RMSF.
For residues whose raw ingredients are available (variants of uncertain
significance at other positions of the same proteins) each protein's
rescaling map is recovered exactly by linear regression and then inverted,
which returns conservation, Norm_Dynamics and mean RMSF for every scored
residue.  Reconstruction is validated on residues where both are known.

Outputs
  cohort_full.csv              labelled cohort, all features and components
  clinvar_all_components.csv   all 768 ClinVar-mapped variants (label-free use)
  protein_scaling_maps.csv     per-protein rescaling parameters
  AdditionalFile2_cohort.csv   public version of the labelled cohort
"""
import json

import numpy as np
import pandas as pd

from depbench.common import (DATA_DIR, PATHOGENIC, BENIGN, add_aa_features,
                             build_cohort, out, read_clinvar_scores)

lab = build_cohort()
allv = add_aa_features(read_clinvar_scores())
print(f"labelled cohort: {len(lab)} variants, {lab.Gene_Name.nunique()} genes, "
      f"{lab.label.sum()} P/LP, {(lab.label == 0).sum()} B/LB")
print(f"all ClinVar-mapped variants: {len(allv)} across {allv.Gene_Name.nunique()} proteins")

# ── recover per-protein rescaling maps from residues with raw ingredients ──
ing = pd.concat([pd.read_csv(f"{DATA_DIR}/VUS_Clinical_Reclassification.csv"),
                 pd.read_csv(f"{DATA_DIR}/case_study_NARS1_PGK1.csv")])
ing = ing.drop_duplicates(["Protein_ID", "Residue_Pos"])
ing["Raw_RCS"] = ing.conservation * (1 - ing.Norm_Dynamics)
assert np.allclose(ing.Raw_DCS, ing.conservation * ing.Norm_Dynamics)
assert np.allclose(np.log1p(ing.DCS_Score) / np.log(11), ing.Normalized_DCS, atol=1e-8)

maps = []
for pid, d in ing.groupby("Protein_ID"):
    d = d.dropna(subset=["Raw_DCS", "Raw_RCS", "DCS_Score", "RCS_Score", "Avg_RMSF"])
    if d.Raw_DCS.nunique() < 2 or d.Raw_RCS.nunique() < 2 or d.Norm_Dynamics.nunique() < 2:
        continue
    a, b = np.polyfit(d.Raw_DCS, d.DCS_Score, 1)
    c, e = np.polyfit(d.Raw_RCS, d.RCS_Score, 1)
    f, h = np.polyfit(d.Norm_Dynamics, d.Avg_RMSF, 1)
    resid = max(np.abs(a * d.Raw_DCS + b - d.DCS_Score).max(),
                np.abs(c * d.Raw_RCS + e - d.RCS_Score).max())
    maps.append(dict(Protein_ID=pid, Gene_Name=d.Gene_Name.iloc[0], n_ref=len(d),
                     dcs_slope=a, dcs_icpt=b, rcs_slope=c, rcs_icpt=e,
                     rmsf_slope=f, rmsf_icpt=h, max_abs_residual=resid,
                     # implied range of the raw quantity across the whole protein
                     raw_dcs_min=-b / a, raw_dcs_max=(10 - b) / a,
                     raw_rcs_min=-e / c, raw_rcs_max=(10 - e) / c))
maps = pd.DataFrame(maps)
maps.to_csv(out("protein_scaling_maps.csv"), index=False)
print(f"\nrescaling maps recovered for {len(maps)} proteins; "
      f"max residual {maps.max_abs_residual.max():.2e} (score units)")


def reconstruct(df):
    df = df.merge(maps[["Protein_ID", "dcs_slope", "dcs_icpt", "rcs_slope", "rcs_icpt",
                        "rmsf_slope", "rmsf_icpt"]], on="Protein_ID", how="left")
    df["Raw_DCS"] = (df.DCS_Score - df.dcs_icpt) / df.dcs_slope
    df["Raw_RCS"] = (df.RCS_Score - df.rcs_icpt) / df.rcs_slope
    df["conservation"] = df.Raw_DCS + df.Raw_RCS
    df["Norm_Dynamics"] = df.Raw_DCS / df.conservation
    df["Avg_RMSF"] = df.rmsf_slope * df.Norm_Dynamics + df.rmsf_icpt
    df["dcs_scale_factor"] = df.dcs_slope       # per-protein multiplier, score units per raw unit
    return df.drop(columns=["dcs_slope", "dcs_icpt", "rcs_slope", "rcs_icpt",
                            "rmsf_slope", "rmsf_icpt"])


lab = reconstruct(lab)
allv = reconstruct(allv)

# validation against residues where the ingredients are known directly
chk = lab.merge(ing[["Protein_ID", "Residue_Pos", "conservation", "Norm_Dynamics", "Avg_RMSF"]],
                on=["Protein_ID", "Residue_Pos"], suffixes=("", "_true"))
val = dict(n_validated=len(chk),
           max_err_conservation=float((chk.conservation - chk.conservation_true).abs().max()),
           max_err_norm_dynamics=float((chk.Norm_Dynamics - chk.Norm_Dynamics_true).abs().max()),
           max_err_rmsf=float((chk.Avg_RMSF - chk.Avg_RMSF_true).abs().max()),
           n_reconstructed=int(lab.conservation.notna().sum()))
print("reconstruction check on labelled residues with known ingredients:", val)
assert lab.conservation.notna().all(), "reconstruction failed for some labelled variants"
json.dump(val, open(out("reconstruction_validation.json"), "w"), indent=1)

# ── label-free within-protein encodings ──────────────────────────────────
# Reference set: every distinct ClinVar-observed residue of the protein
# (any classification, labels never used).  Position-level scores are
# constant across substitutions at a residue, so each residue counts once.
POS_FEATS = ["DCS_Score", "RCS_Score", "Normalized_DCS", "conservation",
             "Norm_Dynamics", "Avg_RMSF", "AM_Mean_Residue_Score"]
ref = allv.drop_duplicates(["Protein_ID", "Residue_Pos"])


def encode(df, feats, reference, key="Protein_ID", tag="ref"):
    """Within-protein percentile, z-score, centring and min-max relative to a
    reference set of residues of the same protein."""
    df = df.copy()
    for f in feats:
        for pid, idx in df.groupby(key).groups.items():
            r = reference.loc[reference[key] == pid, f].dropna().values
            x = df.loc[idx, f].values
            if len(r) == 0:
                continue
            # mid-rank percentile of x within the reference distribution
            less = (r[None, :] < x[:, None]).sum(1)
            eq = (r[None, :] == x[:, None]).sum(1)
            df.loc[idx, f"{f}__pct_{tag}"] = (less + 0.5 * eq) / len(r)
            sd = r.std(ddof=0)
            df.loc[idx, f"{f}__z_{tag}"] = (x - r.mean()) / sd if sd > 0 else 0.0
            df.loc[idx, f"{f}__ctr_{tag}"] = x - r.mean()
            rng_ = r.max() - r.min()
            df.loc[idx, f"{f}__mm_{tag}"] = 10 * (x - r.min()) / rng_ if rng_ > 0 else 5.0
    return df


lab = encode(lab, POS_FEATS, ref, tag="ref")
allv = encode(allv, POS_FEATS, ref, tag="ref")
# AlphaMissense variant-specific score: reference is every ClinVar-observed
# substitution of the protein (substitution-level score, so not deduplicated)
lab = encode(lab, ["AM_Specific_Score"], allv, tag="ref")
allv = encode(allv, ["AM_Specific_Score"], allv, tag="ref")

# Cohort-internal encodings (reference = labelled variants of the same gene)
coh = lab.copy()
lab = encode(lab, POS_FEATS, coh, tag="coh")

lab["n_in_gene"] = lab.groupby("Gene_Name").label.transform("size")
lab["gene_prevalence"] = lab.groupby("Gene_Name").label.transform("mean")
lab.to_csv(out("cohort_full.csv"), index=False)
allv.to_csv(out("clinvar_all_components.csv"), index=False)

# ── public cohort file (Additional file 2) ───────────────────────────────
pub = pd.DataFrame({
    "gene": lab.Gene_Name, "uniprot_accession": lab.UniProt_ID,
    "atlas_protein_id": lab.Protein_ID, "clinvar_variation_id": lab.ClinVar_Accession,
    "protein_change": lab.Protein_Change, "residue_position": lab.Residue_Pos,
    "ref_aa": lab.ref_aa, "alt_aa": lab.alt_aa,
    "clinvar_classification": lab.ClinVar_Significance, "pathogenic": lab.label,
    "n_variants_in_gene": lab.n_in_gene, "gene_pathogenic_fraction": lab.gene_prevalence.round(4),
    "DCS_Score": lab.DCS_Score, "RCS_Score": lab.RCS_Score, "Normalized_DCS": lab.Normalized_DCS,
    "conservation_reconstructed": lab.conservation.round(6),
    "norm_dynamics_reconstructed": lab.Norm_Dynamics.round(6),
    "mean_rmsf_reconstructed_A": lab.Avg_RMSF.round(4),
    "dcs_protein_scale_factor": lab.dcs_scale_factor.round(4),
    "DCS_within_protein_percentile": lab["DCS_Score__pct_ref"].round(4),
    "RCS_within_protein_percentile": lab["RCS_Score__pct_ref"].round(4),
    "AM_variant_score": lab.AM_Specific_Score, "AM_residue_mean_score": lab.AM_Mean_Residue_Score,
})
for c in ["dHydro", "dCharge", "dSize", "dPolarity", "dArom", "dFlex", "absHydro",
          "absCharge", "absSize", "ref_hydro", "ref_flex", "alt_hydro", "alt_flex",
          "ref_is_hinge", "alt_is_hinge", "hinge_disr"]:
    pub[c] = lab[c].round(4)
pub = pub.sort_values(["gene", "residue_position", "protein_change"]).reset_index(drop=True)
pub.to_csv(out("AdditionalFile2_cohort.csv"), index=False)

# consistency with the cohort file distributed with the earlier analysis
orig_path = f"{DATA_DIR}/AdditionalFile2_cohort_v1.csv"
try:
    o = pd.read_csv(orig_path)
    m = o.merge(lab, left_on=["gene", "protein_change"], right_on=["Gene_Name", "Protein_Change"])
    print(f"\nmatch with distributed cohort file: {len(m)}/{len(o)} rows; "
          f"labels identical: {(m.pathogenic == m.label).all()}; "
          f"max |dDCS| = {(m.DCS_Score_x - m.DCS_Score_y).abs().max():.1e}")
except FileNotFoundError:
    pass

print("\nper-gene composition:")
print(lab.groupby("Gene_Name").label.agg(["size", "sum", "mean"])
      .sort_values("size", ascending=False).to_string())
