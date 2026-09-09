"""
extract_slim_stats.py
---------------------
Single-pass extraction from SLiM VCF replicates.
Produces:
  1. {run}_av_DAF.txt      — mean DAF per population, mutation type, DAF bin
  2. {run}_GC.txt           — GC% per population, DAF bin
  3. {run}_SFS.txt          — SFS (DAF histogram) per population

Binning strategy:
  - av_DAF and GC use GLOBAL DAF (pooled across all populations) to assign
    each site to a bin, ensuring each variant appears in exactly one bin
    regardless of population.
  - SFS uses POPULATION-SPECIFIC DAF, reflecting each population's own
    allele frequency distribution.

Usage (specify SLiM run directory:
    python3 extract_slim_stats.py --run /path/to/slimrun --base /path/to/output --nboot 2000

Populations: p1=AFR, p2=EUR, p3=EAS
"""

import numpy as np
import glob
import re
import argparse
import os
import gzip
import pandas as pd
from collections import defaultdict

# ─────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────
S = {"G", "C"}
W = {"A", "T"}

POP_LABELS = {"p1": "AFR", "p2": "EUR", "p3": "EAS"}

DAF_BINS_DEF = {
    "rare-low":  lambda daf: daf < 0.05,
    "rare-high": lambda daf: daf > 0.95,
    "mid":       lambda daf: 0.05 <= daf <= 0.95,
    "all":       lambda daf: True,
}

# GC bins are identical — both keyed on global DAF
GC_BINS_DEF = {
    "rare-low":  lambda daf: daf < 0.05,
    "rare-high": lambda daf: daf > 0.95,
    "mid":       lambda daf: 0.05 <= daf <= 0.95,
    "all":       lambda daf: True,
}

SFS_BINS = np.linspace(0, 1, 21)   # 20 bins of width 0.05
SFS_BIN_LABELS = [
    f"[{SFS_BINS[i]:.2f},{SFS_BINS[i+1]:.2f}{']}' if i == len(SFS_BINS)-2 else ')'}"
    for i in range(len(SFS_BINS) - 1)
]

# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────
def classify(ref, alt):
    """S>W, W>S, or None (skip S>S / W>W)."""
    if ref in S and alt in W:
        return "S>W"
    elif ref in W and alt in S:
        return "W>S"
    return None


def open_vcf(path):
    """Return a text-mode file handle for plain or gzipped VCF."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def extract_rep_index(filename):
    """
    Extract the replicate integer from a filename like p1_rep1000.vcf.gz.
    Returns the integer, or None if no match.
    """
    m = re.search(r"rep(\d+)", os.path.basename(filename))
    return int(m.group(1)) if m else None


def assign_global_daf_bin(global_daf):
    """Return the name of every bin whose condition the global DAF satisfies."""
    return [name for name, cond in DAF_BINS_DEF.items() if cond(global_daf)]


# ─────────────────────────────────────────────
# Global DAF computation for one replicate
# ─────────────────────────────────────────────
def compute_global_daf(vcf_paths):
    """
    Given a list of per-population VCF paths for the SAME replicate,
    compute the global DAF for every site by pooling derived and total
    allele counts across all populations.

    Returns:
        dict: {(chrom, pos, ref, alt): global_daf}
    """
    # site_key → [sum_derived, sum_total]
    site_counts = defaultdict(lambda: [0, 0])

    for vcf in vcf_paths:
        with open_vcf(vcf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 10:
                    continue

                chrom = fields[0]
                pos   = fields[1]
                ref   = fields[3].upper()
                alt   = fields[4].upper()

                if len(ref) != 1 or len(alt) != 1:
                    continue
                if ref not in "ACGT" or alt not in "ACGT":
                    continue

                key = (chrom, pos, ref, alt)
                derived = 0
                total   = 0

                for g in fields[9:]:
                    gt = g.split(":")[0].replace("|", "/")
                    if "." in gt:
                        continue
                    for a in gt.split("/"):
                        if a == "1":
                            derived += 1
                        total += 1

                site_counts[key][0] += derived
                site_counts[key][1] += total

    global_daf = {
        key: counts[0] / counts[1]
        for key, counts in site_counts.items()
        if counts[1] > 0
    }
    return global_daf


# ─────────────────────────────────────────────
# Per-population parsers (DAF + GC)
# ─────────────────────────────────────────────
def parse_vcf_daf(vcf_file, global_daf):
    """
    Read one population VCF. For each site:
      - population-specific DAF is computed from this VCF's genotypes
      - bin assignment uses the pre-computed global DAF

    Returns:
        rep_daf: {bin_name: {mut_type: [pop_daf, ...]}}
        sfs_dafs: [pop_daf, ...]   (for SFS — population-specific, no global binning)
    """
    rep_daf  = {b: {"S>W": [], "W>S": []} for b in DAF_BINS_DEF}
    sfs_dafs = []

    with open_vcf(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos   = fields[1]
            ref   = fields[3].upper()
            alt   = fields[4].upper()

            if len(ref) != 1 or len(alt) != 1:
                continue
            if ref not in "ACGT" or alt not in "ACGT":
                continue

            key = (chrom, pos, ref, alt)
            if key not in global_daf:
                continue   # site absent from global lookup — should not happen

            g_daf = global_daf[key]

            # population-specific DAF
            derived = 0
            total   = 0
            for g in fields[9:]:
                gt = g.split(":")[0].replace("|", "/")
                if "." in gt:
                    continue
                for a in gt.split("/"):
                    if a == "1":
                        derived += 1
                    total += 1

            if total == 0:
                continue

            pop_daf = derived / total
            sfs_dafs.append(pop_daf)

            mut_type = classify(ref, alt)
            if mut_type is None:
                continue

            # bin by GLOBAL DAF
            for bin_name in assign_global_daf_bin(g_daf):
                rep_daf[bin_name][mut_type].append(pop_daf)

    return rep_daf, sfs_dafs


def parse_vcf_gc(vcf_file, global_daf):
    """
    Read one population VCF. Accumulate allele base counts per DAF bin,
    where bin assignment uses the pre-computed global DAF.

    Returns:
        gc_stats: {bin_name: {'A': int, 'C': int, 'G': int, 'T': int}}
    """
    gc_stats = {b: {"A": 0, "C": 0, "G": 0, "T": 0} for b in GC_BINS_DEF}

    with open_vcf(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos   = fields[1]
            ref   = fields[3].upper()
            alt   = fields[4].upper()

            if len(ref) != 1 or len(alt) != 1:
                continue
            if ref not in "ACGT" or alt not in "ACGT":
                continue

            key = (chrom, pos, ref, alt)
            if key not in global_daf:
                continue

            g_daf = global_daf[key]
            bins_for_site = [name for name, cond in GC_BINS_DEF.items() if cond(g_daf)]

            allele_counts = {"A": 0, "C": 0, "G": 0, "T": 0}
            total = 0

            for g in fields[9:]:
                gt = g.split(":")[0].replace("|", "/")
                if "." in gt:
                    continue
                for a in gt.split("/"):
                    if a == "0":
                        allele_counts[ref] += 1
                        total += 1
                    elif a == "1":
                        allele_counts[alt] += 1
                        total += 1

            if total == 0:
                continue

            for bin_name in bins_for_site:
                for base in "ACGT":
                    gc_stats[bin_name][base] += allele_counts[base]

    return gc_stats


# ─────────────────────────────────────────────
# Bootstrap CI
# ─────────────────────────────────────────────
def bootstrap_ci(data, n_boot=1000, ci=95):
    data = np.array(data, dtype=float)
    data = data[~np.isnan(data)]
    if len(data) == 0:
        return np.nan, np.nan, np.nan
    boot_means = [np.mean(np.random.choice(data, size=len(data), replace=True))
                  for _ in range(n_boot)]
    lo = (100 - ci) / 2
    return np.mean(data), np.percentile(boot_means, lo), np.percentile(boot_means, 100 - lo)


# ─────────────────────────────────────────────
# Group VCFs by replicate index
# ─────────────────────────────────────────────
def group_replicates(base_dir, ext):
    """
    Build a dict: rep_index → {pop_key: filepath}
    Only include replicates where ALL populations are present.
    """
    rep_map = defaultdict(dict)

    for pop_key in POP_LABELS:
        pattern = os.path.join(base_dir, f"{pop_key}_rep*.{ext}")
        for path in sorted(glob.glob(pattern)):
            idx = extract_rep_index(path)
            if idx is not None:
                rep_map[idx][pop_key] = path

    n_pops = len(POP_LABELS)
    complete = {idx: paths for idx, paths in rep_map.items() if len(paths) == n_pops}
    incomplete = [idx for idx, paths in rep_map.items() if len(paths) != n_pops]

    if incomplete:
        print(f"  WARNING: skipping {len(incomplete)} replicate(s) with missing population VCFs: "
              f"{sorted(incomplete)}")

    return complete


# ─────────────────────────────────────────────
# Main processing loop
# ─────────────────────────────────────────────
def run_all(base_dir, ext, n_boot=2000):
    """
    Process all replicates. For each replicate:
      1. Compute global DAF across all populations.
      2. For each population, compute pop-specific DAF (binned by global DAF)
         and GC% (binned by global DAF), plus pop-specific SFS.

    Returns DataFrames: (all_daf, all_gc, all_sfs)
    """
    rep_map = group_replicates(base_dir, ext)
    if not rep_map:
        raise FileNotFoundError(f"No complete replicates found in: {base_dir}")

    print(f"  Found {len(rep_map)} complete replicates.")

    # Accumulators: keyed by pop_key
    daf_acc = {p: {b: {"S>W": [], "W>S": []} for b in DAF_BINS_DEF} for p in POP_LABELS}
    daf_n   = {p: {b: {"S>W": [], "W>S": []} for b in DAF_BINS_DEF} for p in POP_LABELS}
    gc_acc  = {p: {b: [] for b in GC_BINS_DEF} for p in POP_LABELS}
    sfs_all = {p: [] for p in POP_LABELS}

    for rep_idx in sorted(rep_map):
        vcf_by_pop = rep_map[rep_idx]

        # ── Step 1: global DAF for this replicate ──
        all_vcfs = list(vcf_by_pop.values())
        g_daf = compute_global_daf(all_vcfs)

        # ── Step 2: per-population stats ──
        for pop_key, vcf_path in vcf_by_pop.items():

            # DAF + SFS
            rep_daf, sfs_dafs = parse_vcf_daf(vcf_path, g_daf)
            sfs_all[pop_key].extend(sfs_dafs)

            for b in DAF_BINS_DEF:
                for k in ["S>W", "W>S"]:
                    vals = rep_daf[b][k]
                    daf_acc[pop_key][b][k].append(np.mean(vals) if vals else np.nan)
                    daf_n[pop_key][b][k].append(len(vals))

            # GC
            gc = parse_vcf_gc(vcf_path, g_daf)
            for bin_name in GC_BINS_DEF:
                c = gc[bin_name]
                total = sum(c.values())
                gc_val = (c["G"] + c["C"]) / total if total > 0 else np.nan
                gc_acc[pop_key][bin_name].append(gc_val)

    # ── Summarise ──
    daf_rows, gc_rows, sfs_rows = [], [], []

    for pop_key, pop_label in POP_LABELS.items():

        # DAF
        for b in DAF_BINS_DEF:
            for k in ["S>W", "W>S"]:
                mean, ci_lo, ci_hi = bootstrap_ci(daf_acc[pop_key][b][k], n_boot)
                daf_rows.append({
                    "Population": pop_label,
                    "Bin":        b,
                    "Mutation":   k,
                    "Mean_DAF":   mean,
                    "CI_lower":   ci_lo,
                    "CI_upper":   ci_hi,
                    "Total_SNPs": int(np.nansum(daf_n[pop_key][b][k])),
                    "N_reps":     len(daf_acc[pop_key][b][k]),
                })

        # GC
        for b in GC_BINS_DEF:
            mean, ci_lo, ci_hi = bootstrap_ci(gc_acc[pop_key][b], n_boot)
            gc_rows.append({
                "Population": pop_label,
                "DAF_bin":    b,
                "Mean_GC":    mean,
                "CI_lower":   ci_lo,
                "CI_upper":   ci_hi,
                "N_reps":     len(gc_acc[pop_key][b]),
            })

        # SFS (population-specific DAF)
        counts, _ = np.histogram(sfs_all[pop_key], bins=SFS_BINS)
        total = counts.sum()
        density = counts / total if total > 0 else counts
        for i, label in enumerate(SFS_BIN_LABELS):
            sfs_rows.append({
                "Population": pop_label,
                "DAF_bin":    label,
                "Count":      int(counts[i]),
                "Proportion": density[i],
            })

    return (
        pd.DataFrame(daf_rows),
        pd.DataFrame(gc_rows),
        pd.DataFrame(sfs_rows),
    )


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Extract DAF, GC%, and SFS from SLiM VCF replicates. "
                    "av_DAF and GC use global DAF for binning; SFS uses population-specific DAF."
    )
    parser.add_argument("--run",   required=True, help="Run name, e.g. run_9")
    parser.add_argument("--base",  required=True, help="Base directory containing VCFs")
    parser.add_argument("--nboot", type=int, default=1000, help="Bootstrap iterations (default: 1000)")
    parser.add_argument("--gz",    action="store_true", help="VCFs are gzipped (.vcf.gz)")
    args = parser.parse_args()

    ext = "vcf.gz" if args.gz else "vcf"
    outdir = args.base
    os.makedirs(outdir, exist_ok=True)

    print(f"\nProcessing run: {args.run}")
    print(f"Base directory: {args.base}")
    print(f"Binning strategy: global DAF for av_DAF + GC; population-specific DAF for SFS\n")

    all_daf, all_gc, all_sfs = run_all(args.base, ext, n_boot=args.nboot)

    daf_out = os.path.join(outdir, f"{args.run}_av_DAF.txt")
    gc_out  = os.path.join(outdir, f"{args.run}_GC.txt")
    sfs_out = os.path.join(outdir, f"{args.run}_SFS.txt")

    all_daf.to_csv(daf_out, sep="\t", index=False)
    all_gc.to_csv(gc_out,   sep="\t", index=False)
    all_sfs.to_csv(sfs_out, sep="\t", index=False)

    print(f"\n✓ Written:")
    print(f"  {daf_out}")
    print(f"  {gc_out}")
    print(f"  {sfs_out}")

    print("\n--- DAF table preview ---")
    print(all_daf.to_string(index=False))
    print("\n--- GC table preview ---")
    print(all_gc.to_string(index=False))


if __name__ == "__main__":
    main()

