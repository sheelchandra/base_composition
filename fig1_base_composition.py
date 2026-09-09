"""
Figure 1: base composition (allele proportions) across populations.

Fig1A - [A] vs [C] proportions across polymorphic sites (common + rare
        variants), colored by AFR vs non-AFR.
Fig1B - GC content per population across fine-grained MAF bins, point size
        scaled by the number of SNPs in each bin.
Fig1C - GC content per population across MAF thresholds (no filter,
        MAF>=1/3/5/10%).

Run with: python fig1_base_composition.py (from inside scripts/)
"""
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns

INPUT_DIR = "/Users/sheelchandra/Downloads/base_comp_github/base_comp_input_files/"
FIG_DIR = "/Users/sheelchandra/Downloads/base_comp_github/base_comp_figs/main_figs/"

# ===========================================================================
# Fig1A
# ===========================================================================
df_all = pd.read_csv(INPUT_DIR + "nucleotide_proportions_all.txt", sep="\t")
df_common = pd.read_csv(INPUT_DIR + "nucleotide_proportions_common.txt", sep="\t")
df_rare = pd.read_csv(INPUT_DIR + "nucleotide_proportions_rare.txt", sep="\t")

pop = pd.read_csv(INPUT_DIR + "20130606_g1k_3202_samples_ped_population.txt", sep=r"\s+")
pop = pop.drop(columns=["FamilyID", "FatherID", "MotherID", "Sex"])
pop = pop.rename(columns={"SampleID": "Sample"})

df_common = df_common.merge(pop, left_on="SAMPLE", right_on="Sample", how="left").drop(columns="Sample")
df_all = df_all.merge(pop, left_on="SAMPLE", right_on="Sample", how="left").drop(columns="Sample")
df_rare = df_rare.merge(pop, left_on="SAMPLE", right_on="Sample", how="left").drop(columns="Sample")

for df in (df_common, df_all, df_rare):
    df["Total"] = df["A"] + df["T"] + df["C"] + df["G"]
    df["prop_A"] = df["A"] / df["Total"]
    df["prop_C"] = df["C"] / df["Total"]

is_afr_common = df_common["Superpopulation"] == "AFR"
is_afr_rare = df_rare["Superpopulation"] == "AFR"

fig, ax = plt.subplots(figsize=(7, 6))

ax.scatter(df_common.loc[~is_afr_common, "prop_A"], df_common.loc[~is_afr_common, "prop_C"],
           color="#3A6EA5", alpha=0.5, s=15, label="Non-AFR")
ax.scatter(df_common.loc[is_afr_common, "prop_A"], df_common.loc[is_afr_common, "prop_C"],
           color="#F73939", alpha=0.5, s=15, label="AFR")
ax.scatter(df_rare.loc[~is_afr_rare, "prop_A"], df_rare.loc[~is_afr_rare, "prop_C"],
           color="#3A6EA5", alpha=0.5, s=15, label="Non-AFR")
ax.scatter(df_rare.loc[is_afr_rare, "prop_A"], df_rare.loc[is_afr_rare, "prop_C"],
           color="#F73939", alpha=0.5, s=15, label="AFR")

ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
ax.spines[["bottom", "left", "top", "right"]].set_linewidth(2)
ax.tick_params(axis="both", which="major", labelsize=17, width=3)
ax.set_xlabel("[A] across polymorphic sites", fontsize=18)
ax.set_ylabel("[C] across polymorphic sites", fontsize=18)  # typo fixed: "acros" -> "across"

plt.tight_layout()
plt.savefig(FIG_DIR + "Fig1A.png", dpi=300)
plt.close()

# ===========================================================================
# Fig1B
# ===========================================================================
df = pd.read_csv(INPUT_DIR + "avg_gc_by_maf_bin_fine_1KG.txt", sep=r"\s+")

custom_colors = ["#00b4d8", "#03045e", "#0077b6", "#90e0ef", "#F73939"]
custom_palette = sns.color_palette(custom_colors)
population_order = ["EAS", "EUR", "AMR", "SAS", "AFR"]

max_maf = 0.20
df["maf_lower"] = df["MAF_bin"].str.split("-").str[0].astype(float)
df = df[df["maf_lower"] < max_maf].sort_values("maf_lower").reset_index(drop=True)
maf_order = df["MAF_bin"].tolist()

maf_labels = [
    f"{lo * 100:.1f}-{hi * 100:.1f}%"
    for lo, hi in (tuple(map(float, b.split("-"))) for b in maf_order)
]
df = df.set_index("MAF_bin").loc[maf_order]

# Scale marker size by N (sqrt scaling so area ~ N)
all_n = np.concatenate([df[f"N_{pop}"].values for pop in population_order])
n_min, n_max = all_n.min(), all_n.max()
size_min, size_max = 20, 600


def size_from_n(n_values):
    frac = (np.sqrt(n_values) - np.sqrt(n_min)) / (np.sqrt(n_max) - np.sqrt(n_min))
    frac = np.clip(frac, 0, 1)
    return size_min + frac * (size_max - size_min)


fig, ax = plt.subplots(figsize=(14, 5.5))
for i, pop in enumerate(population_order):
    x = np.arange(len(maf_order))
    y = df[f"GC_{pop}"].values
    n = df[f"N_{pop}"].values
    ax.plot(x, y, color=custom_palette[i], linewidth=2, zorder=2, label=pop)
    ax.scatter(x, y, s=size_from_n(n), color=custom_palette[i], edgecolor="white", linewidth=0.5, zorder=3)

ax.set_xticks(range(len(maf_order)))
ax.set_xticklabels(maf_labels, fontsize=9, rotation=45, ha="right")
ax.set_xlabel("MAF Bin", fontsize=18)
ax.set_ylabel("GC Content per Population", fontsize=18)
ax.tick_params(axis="y", labelsize=12)
ax.spines[["bottom", "left", "top", "right"]].set_linewidth(2)
ax.tick_params(axis="both", which="major", labelsize=16, width=3)

color_legend = ax.legend(fontsize=13, frameon=False, loc="upper right")
ax.add_artist(color_legend)


def round_down_nice(x):
    """Round x down to the nearest 1/2/5 x 10^k (e.g. 8,342 -> 5,000)."""
    if x <= 0:
        return 0
    exponent = np.floor(np.log10(x))
    magnitude = 10 ** exponent
    for m in (5, 2, 1):
        if x >= m * magnitude:
            return int(m * magnitude)
    return int(magnitude)


example_n_raw = np.geomspace(max(n_min, 1), n_max, num=3)
example_n = sorted(set(round_down_nice(v) for v in example_n_raw))

size_handles = [
    plt.scatter([], [], s=size_from_n(np.array([n]))[0], color="gray", edgecolor="white", label=f"N > {n:,}")
    for n in example_n
]
ax.legend(handles=size_handles, fontsize=11, frameon=False, loc="upper center",
          title="Number of SNPs", title_fontsize=12, labelspacing=1.5, borderpad=1.2)

plt.tight_layout()
plt.savefig(FIG_DIR + "Fig1B.png", dpi=300)
plt.close()

# ===========================================================================
# Fig1C
# ===========================================================================
df = pd.read_csv(INPUT_DIR + "gc_by_maf_1KG.txt", sep="\t", header=None)
df.columns = ["MAF_bin", "AFR", "EUR", "AMR", "EAS", "SAS"]
custom_colors = ["#00b4d8", "#03045e", "#0077b6", "#90e0ef", "#F73939"]
custom_palette = sns.color_palette(custom_colors)
population_order = ["EAS", "EUR", "AMR", "SAS", "AFR"]

maf_order = ["maf0", "maf1", "maf3", "maf5", "maf10"]
maf_labels = ["None", "MAF≥1%", "MAF≥3%", "MAF≥5%", "MAF≥10%"]

df = df.set_index("MAF_bin").loc[maf_order]

fig, ax = plt.subplots(figsize=(7, 6))
for i, pop in enumerate(population_order):
    ax.plot(range(len(maf_order)), df[pop].values, color=custom_palette[i],
            marker="o", label=pop, linewidth=2, markersize=6)

ax.set_xticks(range(len(maf_order)))
ax.set_xticklabels(maf_labels, fontsize=12)
ax.set_xlabel("MAF Threshold", fontsize=18)
ax.set_ylabel("GC Content per Population", fontsize=18)
ax.tick_params(axis="y", labelsize=12)
ax.legend(fontsize=13, frameon=False)
ax.spines[["bottom", "left", "top", "right"]].set_linewidth(2)
ax.tick_params(axis="both", which="major", labelsize=17, width=3)

plt.tight_layout()
plt.savefig(FIG_DIR + "Fig1C.png", dpi=300)
plt.close()

print(f"Wrote Fig1A.png, Fig1B.png, Fig1C.png to {FIG_DIR}")
