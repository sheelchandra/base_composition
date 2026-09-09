"""
Figure 4: comparing empirical GC content to gBGC simulations.

Fig4A - GC content per population, faceted by frequency stratum
        (Rare/Common/All), across three data sources: 1000 Genomes,
        mutation-bias+gBGC simulation, mutation-bias-only (no gBGC)
        simulation.
Fig4B - GC content per population across MAF thresholds, same three
        data sources.

Run with: python fig4_gbgc_sim.py (from inside scripts/)
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

INPUT_DIR = "/Users/sheelchandra/Downloads/base_comp_github/base_comp_input_files/"
FIG_DIR = "/Users/sheelchandra/Downloads/base_comp_github/base_comp_figs/main_figs/"

custom_colors = {
    "AFR": "#F73939",
    "EAS": "#00b4d8",
    "EUR": "#03045e",
    "AMR": "#0077b6",
    "SAS": "#90e0ef",
}


def weighted_gc_by_bin(df):
    """
    Weighted-average GC% (and CI) per Population x Bin, pooling S>W and
    W>S contributions weighted by SNP count, then combining the
    'rare-low' and 'rare-high' bins into a single 'rare' bin.

    Used for both the gBGC-simulation and no-gBGC-control inputs below
    (one shared implementation instead of the original notebook's two
    near-identical copies of this block).
    """
    results = []
    for pop, pop_df in df.groupby("Population"):
        for bin_label, bin_df in pop_df.groupby("Bin"):
            sw = bin_df[bin_df["Mutation"] == "S>W"].iloc[0]
            ws = bin_df[bin_df["Mutation"] == "W>S"].iloc[0]

            # GC contribution: S>W -> DAF, W>S -> 1 - DAF
            gc_sw = 1 - sw["Mean_DAF"]
            gc_ws = ws["Mean_DAF"]

            n_sw = sw["Total_SNPs"]
            n_ws = ws["Total_SNPs"]

            mean_gc = ((gc_sw * n_sw) + (gc_ws * n_ws)) / (n_sw + n_ws)
            ci_low = (ws["CI_lower"] * n_ws + (1 - sw["CI_upper"]) * n_sw) / (n_sw + n_ws)
            ci_high = (ws["CI_upper"] * n_ws + (1 - sw["CI_lower"]) * n_sw) / (n_sw + n_ws)

            results.append({
                "Population": pop,
                "Bin": bin_label,
                "mean_GC": mean_gc,
                "CI_low": ci_low,
                "CI_high": ci_high,
                "n_SNPs": n_sw + n_ws,
            })

    results_df = pd.DataFrame(results)

    rare = results_df[results_df["Bin"].isin(["rare-low", "rare-high"])]
    non_rare = results_df[~results_df["Bin"].isin(["rare-low", "rare-high"])]

    rare_combined = (
        rare.groupby("Population")
        .apply(lambda g: pd.Series({
            "Bin": "rare",
            "mean_GC": np.average(g["mean_GC"], weights=g["n_SNPs"]),
            "CI_low": np.average(g["CI_low"], weights=g["n_SNPs"]),
            "CI_high": np.average(g["CI_high"], weights=g["n_SNPs"]),
            "n_SNPs": g["n_SNPs"].sum(),
        }), include_groups=False)
        .reset_index()
    )

    combined = pd.concat([non_rare, rare_combined], ignore_index=True)
    combined["Bin"] = pd.Categorical(combined["Bin"], categories=["rare", "mid", "all"], ordered=True)
    combined = combined.sort_values(["Population", "Bin"]).reset_index(drop=True)
    return combined


# ===========================================================================
# Fig4A
# ===========================================================================
av_GC_1KG = pd.read_csv(INPUT_DIR + "bootstrap_av_GC.txt", sep="\t")

df = pd.read_csv(INPUT_DIR + "slim_control_av_DAF.txt", sep="\t")
av_GC_slim_no_gBGC = weighted_gc_by_bin(df)

df = pd.read_csv(INPUT_DIR + "slim_gBGC_av_DAF.txt", sep="\t")
av_GC_slim_gBGC = weighted_gc_by_bin(df)

# --- Standardize bin/stratum labels across all three dataframes ---
bin_label_map = {
    "rare": "Rare",
    "mid": "Common",
    "all": "All",
    "DAF < 5% & DAF>95%": "Rare",
    "5%<=DAF<=95%": "Common",
    "All": "All",
}
stratum_order = ["Rare", "Common", "All"]


def prep_sim(df):
    df = df.copy()
    df["Stratum"] = df["Bin"].map(bin_label_map)
    return df[["Population", "Stratum", "mean_GC", "CI_low", "CI_high"]]


def prep_real(df):
    df = df.copy()
    df["Stratum"] = df["Stratum"].map(bin_label_map)
    return df[["Population", "Stratum", "mean_GC", "CI_low", "CI_high"]]


datasets = [
    ("1000 Genomes", prep_real(av_GC_1KG), ["AFR", "AMR", "EUR", "EAS", "SAS"]),
    ("Mutation bias + gBGC", prep_sim(av_GC_slim_gBGC), ["AFR", "EUR", "EAS"]),
    ("Mutation bias, no gBGC", prep_sim(av_GC_slim_no_gBGC), ["AFR", "EUR", "EAS"]),
]
centers = []
ranges = []

for title, df, pops in datasets:
    lo = df["CI_low"].min()
    hi = df["CI_high"].max()
    centers.append((lo + hi) / 2)
    ranges.append(hi - lo)

x_width = max(ranges) * 1.3

for title, df, pops in datasets:
    for _, row in df.iterrows():
        lo = row["mean_GC"] - row["CI_low"]
        hi = row["CI_high"] - row["mean_GC"]
        if lo < 0 or hi < 0:
            print(f"  Pop={row['Population']} Stratum={row['Stratum']} "
                  f"mean={row['mean_GC']:.6f} CI_low={row['CI_low']:.6f} "
                  f"CI_high={row['CI_high']:.6f} xerr_lo={lo:.6f} xerr_hi={hi:.6f}")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax, (title, df, pops) in zip(axes, datasets):
    y_positions = {s: i for i, s in enumerate(stratum_order)}
    n_pops = len(pops)
    offsets = np.linspace(-0.25, 0.25, n_pops)
    center = df["mean_GC"].mean()
    ax.set_xlim(center - x_width / 2, center + x_width / 2)
    for j, pop in enumerate(pops):
        pop_df = df[df["Population"] == pop].set_index("Stratum")

        for stratum in stratum_order:
            if stratum not in pop_df.index:
                continue

            row = pop_df.loc[stratum]
            y = y_positions[stratum] + offsets[j]
            xerr_lo = row["mean_GC"] - row["CI_low"]
            xerr_hi = row["CI_high"] - row["mean_GC"]

            ax.errorbar(row["mean_GC"], y, xerr=[[xerr_lo], [xerr_hi]], fmt="o",
                        color=custom_colors[pop], markersize=7, capsize=3, linewidth=1.2,
                        label=pop if stratum == stratum_order[0] else None)

    ax.set_yticks(range(len(stratum_order)))
    ax.set_yticklabels(stratum_order, fontsize=17)
    ax.set_xlabel("GC content per population", fontsize=20)
    ax.set_title(title, fontsize=22)
    ax.tick_params(axis="x", labelsize=11)
    ax.spines[["bottom", "left", "top", "right"]].set_linewidth(3)
    ax.tick_params(axis="both", which="major", labelsize=17, width=3)

for ax, (title, df, pops), center in zip(axes, datasets, centers):
    ax.set_xlim(center - x_width / 2, center + x_width / 2)
for ax in axes[1:]:
    ax.tick_params(axis="y", labelleft=False)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, fontsize=13, frameon=False, loc="lower center",
           bbox_to_anchor=(0.5, -0.08), ncol=5)

plt.tight_layout()
plt.savefig(FIG_DIR + "Fig4A.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"Wrote Fig4A.png to {FIG_DIR}")

# ===========================================================================
# Fig4B
# ===========================================================================
maf_order = ["maf0", "maf1", "maf3", "maf5", "maf10"]
maf_labels = ["None", "MAF≥1%", "MAF≥3%", "MAF≥5%", "MAF≥10%"]

sim_custom_colors = {
    "p1": "#F73939",  # non-bottlenecked
    "p2": "#03045e",  # bottlenecked
    "p3": "#00b4d8",  # bottlenecked
}
sim_population_order = ["p1", "p2", "p3"]
sim_pop_labels = {"p1": "AFR", "p2": "EUR", "p3": "EAS"}

sim_no_control = pd.read_csv(
    INPUT_DIR + "avg_gc_by_maf_bin_slim.txt", sep=r"\s+"
).set_index("MAF_bin").loc[maf_order]

sim_control = pd.read_csv(
    INPUT_DIR + "avg_gc_by_maf_slim_control_nogBGC.txt", sep=r"\s+"
).set_index("MAF_bin").loc[maf_order]

# Empirical panel: 5 populations, no header row in the file.
emp_custom_colors = ["#00b4d8", "#03045e", "#0077b6", "#90e0ef", "#F73939"]
emp_palette = sns.color_palette(emp_custom_colors)
emp_population_order = ["EAS", "EUR", "AMR", "SAS", "AFR"]

empirical = pd.read_csv(INPUT_DIR + "gc_by_maf_1KG.txt", sep="\t", header=None)
empirical.columns = ["MAF_bin", "AFR", "EUR", "AMR", "EAS", "SAS"]
empirical = empirical.set_index("MAF_bin").loc[maf_order]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
# titles omitted here since Fig4A carries them for the combined figure

ax = axes[1]  # simulation, no control (mutation bias + gBGC)
for pop in sim_population_order:
    ax.plot(range(len(maf_order)), sim_no_control[f"GC_{pop}"].values,
            color=sim_custom_colors[pop], marker="o", label=sim_pop_labels[pop],
            linewidth=2, markersize=6)

ax = axes[2]  # simulation, control (no gBGC)
for pop in sim_population_order:
    ax.plot(range(len(maf_order)), sim_control[f"GC_{pop}"].values,
            color=sim_custom_colors[pop], marker="o", label=sim_pop_labels[pop],
            linewidth=2, markersize=6)

ax = axes[0]  # empirical (1000 Genomes)
for i, pop in enumerate(emp_population_order):
    ax.plot(range(len(maf_order)), empirical[pop].values, color=emp_palette[i],
            marker="o", label=pop, linewidth=2, markersize=6)

for ax in axes:
    ax.set_xticks(range(len(maf_order)))
    ax.set_xticklabels(maf_labels, fontsize=12)
    ax.set_xlabel("MAF Threshold", fontsize=18)
    ax.spines[["bottom", "left", "top", "right"]].set_linewidth(2)
    ax.tick_params(axis="both", which="major", labelsize=14, width=3)
    ax.legend(fontsize=13, frameon=False)

axes[0].set_ylabel("GC Content per Population", fontsize=18)
axes[0].tick_params(axis="y", labelsize=15)

plt.tight_layout()
plt.savefig(FIG_DIR + "Fig4B.png", dpi=300)
plt.close()

print(f"Wrote Fig4B.png to {FIG_DIR}")
