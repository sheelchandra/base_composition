"""
Figure 3: derived allele frequency (DAF) spectrum by mutation type.

Fig3A - broken-y-axis grid of SNP density vs. DAF bin, one panel per
        mutation type (folded across strand-complementary pairs, e.g.
        T>C combined with A>G), colored by population.
Fig3B - mean DAF per mutation type and population, faceted by frequency
        stratum (rare / mid / extreme / all), from bootstrap CIs.

Run with: python fig3_sfs.py (from inside scripts/)
"""
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

INPUT_DIR = "/Users/sheelchandra/Downloads/base_comp_github/base_comp_input_files/"
FIG_DIR = "/Users/sheelchandra/Downloads/base_comp_github/base_comp_figs/main_figs/"

# ===========================================================================
# dataset for population-specific DAF (incl. singletons)
# ===========================================================================
# Column order follows the file's own header (AFR, EUR, AMR, EAS, SAS).
# The original notebook mislabeled this -- it swapped the EUR and AMR
# column names -- which has been fixed here.
#
# dtype="category" for MUTATION avoids materializing ~36M individual Python
# string objects (there are only 12 distinct mutation types) -- this is the
# difference between this ~1.6GB, 36M-row file loading in a couple GB of
# memory versus not fitting in memory at all on a modest machine.
masterdf = pd.read_csv(
    INPUT_DIR + "all_chrs_polarized.txt.gz", sep="\t", header=0,
    dtype={"MUTATION": "category", "DAF_AFR": "float32", "DAF_EUR": "float32",
           "DAF_AMR": "float32", "DAF_EAS": "float32", "DAF_SAS": "float32"},
)  # pandas reads gzip transparently (compression is inferred from ".gz")
masterdf["MUTATION"] = masterdf["MUTATION"].str.upper().astype("category")
masterdf.columns = ["Mutation", "DAF_AFR", "DAF_EUR", "DAF_AMR", "DAF_EAS", "DAF_SAS"]

T_C = masterdf[(masterdf["Mutation"] == "T>C") | (masterdf["Mutation"] == "A>G")]
T_A = masterdf[(masterdf["Mutation"] == "T>A") | (masterdf["Mutation"] == "A>T")]
T_G = masterdf[(masterdf["Mutation"] == "T>G") | (masterdf["Mutation"] == "A>C")]
C_G = masterdf[(masterdf["Mutation"] == "C>G") | (masterdf["Mutation"] == "C>G")]
C_T = masterdf[(masterdf["Mutation"] == "C>T") | (masterdf["Mutation"] == "G>A")]
C_A = masterdf[(masterdf["Mutation"] == "C>A") | (masterdf["Mutation"] == "G>T")]
W_S = pd.concat([T_C, T_G], ignore_index=True)
W_W = T_A.copy()
S_S = C_G.copy()
S_W = pd.concat([C_T, C_A], ignore_index=True)

# ===========================================================================
# Fig3A
# ===========================================================================
def _compute_density(df, population_order, frequency_bins):
    n_bins = len(frequency_bins) - 1
    total_rows = len(df)
    population_density = {pop: [0.0] * n_bins for pop in population_order}
    for pop in population_order:
        col = f"DAF_{pop}"
        for i in range(n_bins):
            lo, hi = frequency_bins[i], frequency_bins[i + 1]
            if i == n_bins - 1:  # last bin: inclusive of upper limit
                count = len(df[(df[col] >= lo) & (df[col] <= hi)])
            else:
                count = len(df[(df[col] >= lo) & (df[col] < hi)])
            population_density[pop][i] = count / total_rows
    return population_density


def _draw_break_marks(ax_top, ax_bottom):
    d = .015
    kwargs = dict(transform=ax_top.transAxes, color="k", clip_on=False, linewidth=1.5)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    kwargs.update(transform=ax_bottom.transAxes)
    ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)


def plot_sfs_broken_grid(
    dataframes,
    mutation_order=("C>T", "C>A", "C>G", "T>C", "T>G", "T>A"),
    n_rows=3,
    n_cols=2,
    population_order=("EAS", "EUR", "AMR", "SAS", "AFR"),
    bin_width=0.05,
    bottom_ylim=0.03,
    top_ylim_start=0.85,
    custom_palette=None,
    figsize=(15, 15),
    height_ratios=(1, 2),
    save_path=None,
):
    """Grid of broken-y-axis SFS plots, one panel per mutation type."""
    for m in mutation_order:
        if m not in dataframes:
            raise ValueError(f"'{m}' not found in `dataframes`.")
        for pop in population_order:
            col = f"DAF_{pop}"
            if col not in dataframes[m].columns:
                raise ValueError(f"Column '{col}' not found in dataframe for '{m}'.")

    frequency_bins = np.clip(np.arange(0.0, 1.0 + bin_width, bin_width), 0.0, 1.0)
    n_bins = len(frequency_bins) - 1

    if custom_palette is None:
        custom_palette = sns.color_palette("tab10", len(population_order))

    bar_width = 0.8 / len(population_order)
    x = np.arange(n_bins)
    offset = np.linspace(-0.4 + bar_width / 2, 0.4 - bar_width / 2, len(population_order))
    custom_tick_labels = [
        f'[{frequency_bins[i]:.2f}-{frequency_bins[i+1]:.2f}{")" if i < n_bins - 1 else "]"}'
        for i in range(n_bins)
    ]

    fig = plt.figure(figsize=figsize)
    outer = gridspec.GridSpec(n_rows, n_cols, figure=fig, hspace=0.15, wspace=0.12)

    legend_handles, legend_labels = None, None

    for idx, mutation in enumerate(mutation_order):
        row, col = divmod(idx, n_cols)
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[row, col], height_ratios=height_ratios, hspace=0.05
        )
        ax_top = fig.add_subplot(inner[0])
        ax_bottom = fig.add_subplot(inner[1], sharex=ax_top)

        density = _compute_density(dataframes[mutation], population_order, frequency_bins)
        all_values = np.array([density[pop] for pop in population_order])
        overall_max = all_values.max()

        bottom_cut = bottom_ylim.get(mutation, 0.03) if isinstance(bottom_ylim, dict) else bottom_ylim
        top_cut = top_ylim_start.get(mutation, 0.75) if isinstance(top_ylim_start, dict) else top_ylim_start

        for i, pop in enumerate(population_order):
            for ax in (ax_top, ax_bottom):
                ax.bar(x + offset[i], density[pop], width=bar_width, label=pop, color=custom_palette[i])
        if idx == 0:
            legend_handles, legend_labels = ax_bottom.get_legend_handles_labels()

        ax_top.set_ylim(top_cut, overall_max * 1.05)
        ax_bottom.set_ylim(0, bottom_cut)

        ax_top.spines["bottom"].set_visible(False)
        ax_bottom.spines["top"].set_visible(False)
        ax_top.tick_params(bottom=False)
        _draw_break_marks(ax_top, ax_bottom)

        ax_top.set_title(mutation, fontsize=16)

        if row == n_rows - 1:
            ax_bottom.set_xticks(x)
            ax_bottom.set_xticklabels(custom_tick_labels, rotation=45, ha="right", fontsize=9)
        else:
            ax_bottom.set_xticks(x)
            ax_bottom.set_xticklabels([])

        for ax in (ax_top, ax_bottom):
            ax.spines[["left", "right"]].set_linewidth(1.5)
            ax.tick_params(axis="both", which="major", labelsize=10, width=1.2)
        ax_top.spines["top"].set_linewidth(1.5)
        ax_bottom.spines["bottom"].set_linewidth(1.5)

    fig.text(0.05, 0.5, "SNP Density", va="center", rotation="vertical", fontsize=18)
    fig.text(0.5, 0.045, "DAF Bins", ha="center", fontsize=18)

    fig.legend(legend_handles, legend_labels, loc="upper center",
               bbox_to_anchor=(0.5, 0.94), ncol=len(population_order), fontsize=13, frameon=False)

    plt.tight_layout(rect=[0.06, 0.07, 1, 0.93])
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


custom_colors = ["#00b4d8", "#03045e", "#0077b6", "#90e0ef", "#F73939"]
custom_palette = sns.color_palette(custom_colors)

dataframes = {"C>T": C_T, "C>A": C_A, "C>G": C_G, "T>C": T_C, "T>G": T_G, "T>A": T_A}
plot_sfs_broken_grid(
    dataframes,
    custom_palette=custom_palette,
    save_path=FIG_DIR + "Fig3A.png",
)
print(f"Wrote Fig3A.png to {FIG_DIR}")

# ===========================================================================
# Fig3B
# ===========================================================================
results_df = pd.read_csv(INPUT_DIR + "bootstrap_av_DAF.txt", sep="\t", header=0)

population_order = ["AFR", "EAS", "EUR", "AMR", "SAS"]
mut_order = ["C>T", "C>A", "C>G", "T>C", "T>G", "T>A"]
stratum_order = ["DAF<5%", "5%<=DAF<=95%", "DAF>95%", "All"]
stratum_labels = {
    "DAF<5%": "DAF < 5% (N=32,793,679)",
    "5%<=DAF<=95%": "5% ≤ DAF ≤ 95% (N=2,470,918)",
    "DAF>95%": "DAF > 95% (N=527,914)",
    "All": "All SNPs (N=35,792,511)",
}
custom_colors = {
    "AFR": "#F73939",
    "EAS": "#00b4d8",
    "EUR": "#03045e",
    "AMR": "#0077b6",
    "SAS": "#90e0ef",
}

n_pops = len(population_order)
x = np.arange(len(mut_order))
offsets = np.linspace(-0.25, 0.25, n_pops)

fig, axes = plt.subplots(2, 2, figsize=(13, 10))
axes_flat = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]

for idx, (ax, stratum) in enumerate(zip(axes_flat, stratum_order)):
    sub = results_df[results_df["Stratum"] == stratum]

    for j, pop in enumerate(population_order):
        pop_data = sub[sub["Population"] == pop].set_index("MutClass")

        y = [pop_data.loc[m, "mean_DAF"] if m in pop_data.index else np.nan for m in mut_order]
        y_low = [pop_data.loc[m, "CI_low"] if m in pop_data.index else np.nan for m in mut_order]
        y_high = [pop_data.loc[m, "CI_high"] if m in pop_data.index else np.nan for m in mut_order]

        yerr_low = np.array(y) - np.array(y_low)
        yerr_high = np.array(y_high) - np.array(y)

        ax.errorbar(x + offsets[j], y, yerr=[yerr_low, yerr_high], fmt="o",
                    color=custom_colors[pop], markersize=7, capsize=3, linewidth=1.2, label=pop)

    ax.set_xticks(x)
    ax.spines[["bottom", "left", "top", "right"]].set_linewidth(3)
    ax.tick_params(axis="both", which="major", labelsize=17, width=3)

    if idx in (2, 3):
        ax.set_xticklabels(mut_order, fontsize=17)
        ax.set_xlabel("Mutation Type", fontsize=20)
    else:
        ax.set_xticklabels([])

    if idx in (0, 2):
        ax.set_ylabel("Average DAF", fontsize=20)
    else:
        ax.tick_params(axis="y")

    ax.set_title(stratum_labels[stratum], fontsize=22)

handles, labels = axes_flat[0].get_legend_handles_labels()
fig.legend(handles, labels, fontsize=14, frameon=False, loc="lower center",
           bbox_to_anchor=(0.5, -0.04), ncol=n_pops)

plt.tight_layout()
plt.savefig(FIG_DIR + "Fig3B.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"Wrote Fig3B.png to {FIG_DIR}")
