import gzip, sys, collections
from pathlib import Path

# --- Config ---
SAMPLE_INFO   = "non_admixed_1KG_sample_IDs.txt"   # txt file with sample IDs of individuals you want to include
VCF_PATTERN   = "all/chr{chrom}_all.vcf.gz"         # specify variant class, can change to common/ or rare/
CHROMS        = [str(c) for c in range(1, 23)]
CATEGORY      = "all"   # "all", "common", "rare"

VCF_PATTERNS = {
    "all":    "all/chr{chrom}_all.vcf.gz",
    "common": "common/chr{chrom}_common.vcf.gz",
    "rare":   "rare/chr{chrom}_rare.vcf.gz",
}

# --- Load sample metadata ---
sample_pop = {}   # sample_id -> superpopulation (AFR/EUR/AMR/EAS/SAS)
with open(SAMPLE_INFO) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3:
            sample_pop[parts[0]] = parts[2]

# --- Per-individual nucleotide counters ---
# counts[sample_id] = {'A':0,'C':0,'G':0,'T':0}
counts = {s: {'A':0,'C':0,'G':0,'T':0} for s in sample_pop}

def score_site(ref, alt, sample_gts, sample_ids):
    """Add nucleotide contributions for one biallelic SNP site."""
    ref = ref.upper()
    alt = alt.upper()
    if ref not in 'ACGT' or alt not in 'ACGT':
        return  # skip indels or ambiguous

    for sample_id, gt_field in zip(sample_ids, sample_gts):
        if sample_id not in counts:
            continue
        gt = gt_field.split(':')[0]          # grab GT, ignore other FORMAT fields
        alleles = gt.replace('|', '/').split('/')
        if '.' in alleles:
            continue                          # missing genotype
        try:
            n_alt = sum(int(a) for a in alleles)   # 0, 1, or 2
        except ValueError:
            continue
        n_ref = 2 - n_alt
        counts[sample_id][ref] += n_ref
        counts[sample_id][alt] += n_alt

# --- Main loop over chromosomes and categories ---
for category in ['all', 'common', 'rare']:
    # Reset counts per category
    counts = {s: {'A':0,'C':0,'G':0,'T':0} for s in sample_pop}

    for chrom in CHROMS:
        vcf_path = VCF_PATTERNS[category].format(chrom=chrom)
        if not Path(vcf_path).exists():
            print(f"  Skipping missing: {vcf_path}", file=sys.stderr)
            continue

        print(f"  [{category}] chr{chrom}...", file=sys.stderr)

        with gzip.open(vcf_path, 'rt') as f:
            sample_ids = []
            for line in f:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    sample_ids = line.strip().split('\t')[9:]
                    continue
                parts = line.strip().split('\t')
                ref, alt = parts[3], parts[4]
                sample_gts = parts[9:]
                score_site(ref, alt, sample_gts, sample_ids)

    # --- Aggregate by superpopulation ---
    pop_counts = collections.defaultdict(lambda: {'A':0,'C':0,'G':0,'T':0})
    for sample_id, nuc_counts in counts.items():
        pop = sample_pop.get(sample_id)
        if pop:
            for nuc, n in nuc_counts.items():
                pop_counts[pop][nuc] += n

    # --- Write output ---
    outfile = f"nucleotide_proportions_{category}.txt"
    with open(outfile, 'w') as out:
        out.write("category\tpopulation\tA\tC\tG\tT\ttotal\tprop_A\tprop_C\n")
        for pop in sorted(pop_counts):
            c = pop_counts[pop]
            total = sum(c.values())
            if total == 0:
                continue
            out.write(
                f"{category}\t{pop}\t"
                f"{c['A']}\t{c['C']}\t{c['G']}\t{c['T']}\t{total}\t"
                f"{c['A']/total:.6f}\t{c['C']/total:.6f}\n"
            )
    print(f"Written: {outfile}")

