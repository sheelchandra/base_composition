"""

Uses masterdf (polarized, per-population DAF for all 1KG SNPs; produced by
polarize_and_combine.sh) to compute mean DAF for each mutation class x DAF
stratum x population, with block bootstrap confidence intervals.

Output: bootstrap_av_DAF.txt (used for Figure 3B)
"""

import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
INPUT_DIR = "base_comp_input_files/"
OUTDIR    = "base_comp_input_files/"

MASTERDF_FILE = INPUT_DIR + "all_chrs_polarized.txt"
OUT_FILE      = OUTDIR + "bootstrap_av_DAF.txt"

# ---------------------------------------------------------------------------
# Load masterdf: one row per polarized, biallelic SNP with per-population DAF
# ---------------------------------------------------------------------------
masterdf = pd.read_csv(MASTERDF_FILE, sep='\t', header=0)
masterdf['Mutation'] = masterdf['Mutation'].str.upper()
masterdf.columns = ['Mutation', 'DAF_AFR', 'DAF_EUR', 'DAF_AMR', 'DAF_EAS', 'DAF_SAS', 'DAF_global']

# --- Setup ---
# Each key is the label, value is the list of Mutation values it includes
mutation_classes = {
    'C>T': ['C>T', 'G>A'],
    'C>A': ['C>A', 'G>T'],
    'C>G': ['C>G', 'G>C'],
    'T>C': ['T>C', 'A>G'],
    'T>G': ['T>G', 'A>C'],
    'T>A': ['T>A', 'A>T'],
}

populations = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
n_blocks    = 1000
n_resamples = 1000
ci_lower    = 2.5
ci_upper    = 97.5

# --- DAF strata defined on global DAF ---
strata = {
    'DAF<5%'       : masterdf['DAF_global'] < 0.05,
    'DAF>95%'      : masterdf['DAF_global'] > 0.95,
    '5%<=DAF<=95%' : (masterdf['DAF_global'] >= 0.05) & (masterdf['DAF_global'] <= 0.95),
    'All'          : pd.Series(True, index=masterdf.index),
}


def block_bootstrap(values, n_blocks, n_resamples):
    """
    Partition values into n_blocks, resample blocks with replacement
    n_resamples times, return array of replicate means.
    """
    n = len(values)
    block_edges = np.array_split(np.arange(n), n_blocks)  # non-overlapping blocks

    replicate_means = np.empty(n_resamples)
    for r in range(n_resamples):
        sampled_indices = np.concatenate(
            [block_edges[i] for i in np.random.randint(0, n_blocks, size=n_blocks)]
        )
        replicate_means[r] = values[sampled_indices].mean()

    return replicate_means


# --- Run: mean DAF + block-bootstrap CI per mutation class x stratum x population ---
results = []

for mut_label, mut_list in mutation_classes.items():
    mut_mask = masterdf['Mutation'].isin(mut_list)
    print(f"Processing {mut_label} ({mut_mask.sum()} SNPs)...")

    for stratum_label, stratum_mask in strata.items():
        combined = masterdf[mut_mask & stratum_mask]

        for pop in populations:
            col    = f'DAF_{pop}'
            values = combined[col].values

            if len(values) < n_blocks:
                print(f"  Skipping {mut_label} | {stratum_label} | {pop} — too few SNPs ({len(values)})")
                results.append({
                    'MutClass'  : mut_label,
                    'Population': pop,
                    'Stratum'   : stratum_label,
                    'n_SNPs'    : len(values),
                    'mean_DAF'  : np.nan,
                    'CI_low'    : np.nan,
                    'CI_high'   : np.nan,
                })
                continue

            replicates = block_bootstrap(values, n_blocks, n_resamples)

            results.append({
                'MutClass'  : mut_label,
                'Population': pop,
                'Stratum'   : stratum_label,
                'n_SNPs'    : len(values),
                'mean_DAF'  : values.mean(),
                'CI_low'    : np.percentile(replicates, ci_lower),
                'CI_high'   : np.percentile(replicates, ci_upper),
            })

results_df = pd.DataFrame(results)

results_df.to_csv(OUT_FILE, sep='\t', index=False)
print(f"Done: {OUT_FILE}")