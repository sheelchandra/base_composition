"""

Uses masterdf (polarized, per-population DAF for all 1KG SNPs; produced by
polarize_and_combine.sh) to compute GC% (S-derived allele frequency) per
population within DAF strata (rare, common, all variants), with block
bootstrap confidence intervals.

Output: bootstrap_av_GC.txt (used for Figure 4A, 1KG)
"""

import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
INPUT_DIR = "base_comp_input_files/"
OUTDIR    = "base_comp_input_files/"

MASTERDF_FILE = INPUT_DIR + "all_chrs_polarized.txt"
OUT_FILE      = OUTDIR + "bootstrap_av_GC.txt"

# ---------------------------------------------------------------------------
# Load masterdf: one row per polarized, biallelic SNP with per-population DAF
# ---------------------------------------------------------------------------
masterdf = pd.read_csv(MASTERDF_FILE, sep='\t', header=0)
masterdf['Mutation'] = masterdf['Mutation'].str.upper()
masterdf.columns = ['Mutation', 'DAF_AFR', 'DAF_EUR', 'DAF_AMR', 'DAF_EAS', 'DAF_SAS', 'DAF_global']

populations = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
ci_lower    = 2.5
ci_upper    = 97.5

# --- Strata based on global DAF ---
strata = {
    'DAF < 5% & DAF>95%' : (masterdf['DAF_global'] > 0.95) | (masterdf['DAF_global'] < 0.05),
    '5%<=DAF<=95%'       : (masterdf['DAF_global'] >= 0.05) & (masterdf['DAF_global'] <= 0.95),
    'All'                : pd.Series(True, index=masterdf.index),
}

# --- S/W derived mutation classes ---
S_derived = ['C>G', 'T>C', 'T>G', 'A>C', 'A>G', 'G>C']
W_derived = ['C>T', 'C>A', 'T>A', 'G>A', 'G>T', 'A>T']


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


# --- Per-population GC%: DAF for S-derived mutations, 1 - DAF for W-derived ---
gc_df = masterdf.copy()
for pop in populations:
    daf_col = f'DAF_{pop}'
    gc_df[f'GC_{pop}'] = np.where(
        gc_df['Mutation'].isin(S_derived),
        gc_df[daf_col],
        1 - gc_df[daf_col]
    )

# --- Run: mean GC% + block-bootstrap CI per stratum x population ---
results_gc_overall = []

for stratum_label, stratum_mask in strata.items():
    combined = gc_df[stratum_mask]

    for pop in populations:
        col    = f'GC_{pop}'
        values = combined[col].values

        replicates = block_bootstrap(values, 100, 500)

        results_gc_overall.append({
            'Population': pop,
            'Stratum'   : stratum_label,
            'n_SNPs'    : len(values),
            'mean_GC'   : values.mean(),
            'CI_low'    : np.percentile(replicates, ci_lower),
            'CI_high'   : np.percentile(replicates, ci_upper),
        })

results_gc_overall_df = pd.DataFrame(results_gc_overall)

results_gc_overall_df.to_csv(OUT_FILE, sep='\t', index=False)
print(f"Done: {OUT_FILE}")