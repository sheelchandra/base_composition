"""
Figure 2: bottleneck effect on base composition, across four species.

For each species (human, silkworm, maize, mouse), plots [A] vs [C] per
individual, colored by bottleneck status. Each individual contributes two
points -- one from the "Common (MAF>=5%)" SNP set and one from the "All"
SNP set -- so two clusters appear naturally in each panel.

Run with: python fig2_bottleneck.py (from inside scripts/)
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_DIR = "/Users/sheelchandra/Downloads/base_comp_github/base_comp_input_files/"
FIG_DIR = "/Users/sheelchandra/Downloads/base_comp_github/base_comp_figs/main_figs/"


def load_props(path, sep="\t"):
    """Read a *_sums.txt file of per-sample allele counts and convert to proportions."""
    d = pd.read_csv(path, sep=sep, low_memory=False, names=["Sample", "A", "G", "T", "C"], skiprows=1)
    allele_sum = d[["A", "G", "T", "C"]].sum(axis=1)
    props = pd.DataFrame()
    props["Sample"] = d["Sample"]
    for allele in ["A", "G", "T", "C"]:
        props[allele] = d[allele] / allele_sum
    props["prop_GC"] = props["C"] + props["G"]
    return props


# ===========================================================================
# SILKWORM
# ===========================================================================
pop_silk = pd.read_excel(INPUT_DIR + "silkworm sample names.xlsx")

silk_common = load_props(INPUT_DIR + "silkworm_common_sums.txt")
silk_common = silk_common.merge(pop_silk[["Sample", "Population"]], on="Sample", how="inner")

silk_all = load_props(INPUT_DIR + "silkworm_all_sums.txt")
silk_all = silk_all.merge(pop_silk[["Sample", "Population"]], on="Sample", how="inner")

silk_common["Group"] = np.where(silk_common["Population"] == "Wild", "Non-Bottlenecked", "Bottlenecked")
silk_all["Group"] = np.where(silk_all["Population"] == "Wild", "Non-Bottlenecked", "Bottlenecked")

silk_common["VariantSet"] = "Common (MAF≥5%)"
silk_all["VariantSet"] = "All"
silk_df = pd.concat([silk_all, silk_common], ignore_index=True)

# ===========================================================================
# MAIZE
# ===========================================================================
pop_maize = pd.read_excel(INPUT_DIR + "Maize HapMap3 sample names.xlsx")

maize_all = load_props(INPUT_DIR + "maize_all_sums.txt")
maize_all = maize_all.merge(pop_maize, left_on="Sample", right_on="Prefixed_Taxon", how="left")
maize_all.dropna(subset=["Species", "Prefixed_Taxon"], inplace=True)
maize_all.drop(columns=["Prefixed_Taxon"], inplace=True)

maize_common = load_props(INPUT_DIR + "maize_common_sums.txt")
maize_common = maize_common.merge(pop_maize, left_on="Sample", right_on="Prefixed_Taxon", how="left")
maize_common.dropna(subset=["Species", "Prefixed_Taxon"], inplace=True)
maize_common.drop(columns=["Prefixed_Taxon"], inplace=True)

maize_all["Group"] = np.where(maize_all["Species"] == " Z. mays spp. Mays", "Bottlenecked", "Non-Bottlenecked")
maize_common["Group"] = np.where(maize_common["Species"] == " Z. mays spp. Mays", "Bottlenecked", "Non-Bottlenecked")

maize_common["VariantSet"] = "Common (MAF≥5%)"
maize_all["VariantSet"] = "All"
maize_df = pd.concat([maize_all, maize_common], ignore_index=True)

# ===========================================================================
# MOUSE
# ===========================================================================
pop_wild = pd.read_csv(INPUT_DIR + "mouse_samplenames.txt", sep="\t", low_memory=False,
                        names=["Sample"], skiprows=1)
# this analysis does not include spretus individuals, so exclude them
spre_exclude = ["Ms_SPRE1_SP36.variant9", "Ms_SPRE2_SP39.variant9", "Ms_SPRE3_SP41.variant9",
                "Ms_SPRE4_SP51.variant9", "Ms_SPRE5_SP62.variant9", "Ms_SPRE6_SP68.variant9",
                "Ms_SPRE7_SP69.variant9", "Ms_SPRE8_SP70.variant9"]

mouse_common = load_props(INPUT_DIR + "mouse_common_sums.txt")
mouse_common = mouse_common.merge(pop_wild[["Sample"]], on="Sample", how="left")
mouse_common = mouse_common[~mouse_common["Sample"].isin(spre_exclude)]
mouse_common["Population"] = mouse_common["Sample"].str[:3]

mouse_all = load_props(INPUT_DIR + "mouse_all_sums.txt")
mouse_all = mouse_all.merge(pop_wild[["Sample"]], on="Sample", how="left")
mouse_all = mouse_all[~mouse_all["Sample"].isin(spre_exclude)]
mouse_all["Population"] = mouse_all["Sample"].str[:3]

mouse_common["Group"] = np.where(mouse_common["Population"] == "Mmc", "Non-Bottlenecked", "Bottlenecked")
mouse_all["Group"] = np.where(mouse_all["Population"] == "Mmc", "Non-Bottlenecked", "Bottlenecked")

mouse_common["VariantSet"] = "Common (MAF≥5%)"
mouse_all["VariantSet"] = "All"
mouse_df = pd.concat([mouse_all, mouse_common], ignore_index=True)

# ===========================================================================
# HUMAN
# ===========================================================================
pop_human = pd.read_csv(INPUT_DIR + "20130606_g1k_3202_samples_ped_population.txt", sep=r"\s+")
pop_human = pop_human.drop(columns=["FamilyID", "FatherID", "MotherID", "Sex"])
pop_human = pop_human.rename(columns={"SampleID": "Sample"})


def load_human(path):
    d = pd.read_csv(path, sep="\t")
    d["Total"] = d["A"] + d["T"] + d["C"] + d["G"]
    d["prop_GC"] = (d["C"] + d["G"]) / d["Total"]
    d = d.merge(pop_human, left_on="SAMPLE", right_on="Sample", how="left").drop(columns="Sample")
    return d


human_all = load_human(INPUT_DIR + "nucleotide_proportions_all.txt")
human_common = load_human(INPUT_DIR + "nucleotide_proportions_common.txt")

human_all["Group"] = np.where(human_all["Superpopulation"] == "AFR", "Non-Bottlenecked", "Bottlenecked")
human_common["Group"] = np.where(human_common["Superpopulation"] == "AFR", "Non-Bottlenecked", "Bottlenecked")

human_common["VariantSet"] = "Common (MAF≥5%)"
human_all["VariantSet"] = "All"
human_df = pd.concat([human_all, human_common], ignore_index=True)

# ===========================================================================
# PLOT: one subplot per species, arranged 2x2.
# Within each subplot: x = [A], y = [C], one point per individual. Color
# encodes bottleneck status only, not variant set.
# ===========================================================================
species_panels = [
    ("Human (n=2,064)", human_df),
    ("Silkworm (n=1,082)", silk_df),
    ("Maize (n=914)", maize_df),
    ("Mouse (n=59)", mouse_df),
]

group_colors = {
    "Non-Bottlenecked": "#F73939",  # red
    "Bottlenecked": "#3A6EA5",      # blue
}

fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=False, sharey=False)
axes = axes.flatten()

for i, (ax, (title, d)) in enumerate(zip(axes, species_panels)):
    for group, color in group_colors.items():
        sub = d[d["Group"] == group]
        ax.scatter(sub["A"], sub["C"], color=color, alpha=0.75, s=15,
                   label=group if i == 0 else None)
    ax.set_title(title, fontsize=20)
    ax.spines[["bottom", "left", "top", "right"]].set_linewidth(2)
    ax.tick_params(axis="both", which="major", labelsize=16, width=2)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="center", bbox_to_anchor=(0.5, 0.05),
           ncol=2, fontsize=15, markerscale=4, frameon=False)

plt.tight_layout(rect=[0.02, 0.10, 1, 1])
plt.subplots_adjust(hspace=0.2, wspace=0.2)

fig.supxlabel("[A] across polymorphic sites", fontsize=19, y=0.08)
fig.supylabel("[C] across polymorphic sites", fontsize=19, x=0.01)
plt.savefig(FIG_DIR + "Fig2.png", dpi=300)
plt.close()

print(f"Wrote Fig2.png to {FIG_DIR}")
