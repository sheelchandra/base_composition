"""
Figure S3: base composition (allele proportions) in high vs. low
recombination-rate regions, common variants only.

Run with: python figS3_recombination.py (from inside scripts/)
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_DIR = "/Users/sheelchandra/Downloads/base_comp_github/base_comp_input_files/"
FIG_DIR = "/Users/sheelchandra/Downloads/base_comp_github/base_comp_figs/supp_figs/"

df_common_high = pd.read_csv(INPUT_DIR + "nucleotide_proportions_high.txt", sep="\t")
df_common_low = pd.read_csv(INPUT_DIR + "nucleotide_proportions_low.txt", sep="\t")

pop = pd.read_csv(INPUT_DIR + "20130606_g1k_3202_samples_ped_population.txt", sep=r"\s+")
pop = pop.drop(columns=["FamilyID", "FatherID", "MotherID", "Sex"])
pop = pop.rename(columns={"SampleID": "Sample"})

# Merge with population annotations
df_common_low = df_common_low.merge(pop, left_on="SAMPLE", right_on="Sample", how="left").drop(columns="Sample")
df_common_high = df_common_high.merge(pop, left_on="SAMPLE", right_on="Sample", how="left").drop(columns="Sample")


def tag(d, label):
    d = d.copy()
    d["VariantSet"] = label
    return d


df_common_high = tag(df_common_high, "Common (High Recombination)")
df_common_low = tag(df_common_low, "Common (Low Recombination)")

# Combine both variant sets into one long dataframe for plotting
plot_df = pd.concat([df_common_high, df_common_low], ignore_index=True)
plot_df["Group"] = np.where(plot_df["Superpopulation"] == "AFR", "Non-Bottlenecked", "Bottlenecked")

# Plot: x = [A], y = [C], one point per individual per variant set.
# Color encodes Group (bottleneck status).
colors = {"Bottlenecked": "#3A6EA5", "Non-Bottlenecked": "#F73939"}

fig, ax = plt.subplots(figsize=(8, 6))
for group, color in colors.items():
    sub = plot_df[plot_df["Group"] == group]
    ax.scatter(sub["A"], sub["C"], color=color, alpha=0.75, s=15, label=group)

ax.set_xlabel("[A] across polymorphic sites", fontsize=18)
ax.set_ylabel("[C] across polymorphic sites", fontsize=18)
ax.spines[["bottom", "left", "top", "right"]].set_linewidth(2)
ax.tick_params(axis="both", which="major", labelsize=15, width=3)
ax.legend(fontsize=14, markerscale=4, frameon=False)
plt.tight_layout()
plt.savefig(FIG_DIR + "FigS3.png", dpi=300)
plt.close()

print(f"Wrote FigS3.png to {FIG_DIR}")
