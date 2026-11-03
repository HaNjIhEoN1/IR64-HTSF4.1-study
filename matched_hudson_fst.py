#!/usr/bin/env python3
"""
Direct 34,649-bp matched-window Hudson FST analysis.

This script does not use the previous 10-kb sliding-window results.

It:
  1. calculates Hudson FST directly from the VCF for the exact
     qHTSF4.1 candidate LD interval;
  2. divides the genome into non-overlapping 34,649-bp controls;
  3. calculates Hudson FST directly for every control region;
  4. matches controls by the number of informative SNPs;
  5. calculates a genome-wide empirical upper-tail probability.

Candidate interval:
    chr04:18,535,070-18,569,718
    length = 34,649 bp

Empirical P:
    (1 + number of matched controls with FST >= candidate FST)
    / (1 + total number of matched controls)
"""

import argparse
import sys
from pathlib import Path

import allel
import numpy as np
import pandas as pd


# ================================================================
# Fixed coordinates: 1-based inclusive
# ================================================================

CANDIDATE_CHROM = "chr04"
CANDIDATE_START = 18_535_070
CANDIDATE_END = 18_569_718

REGION_LENGTH = (
    CANDIDATE_END
    - CANDIDATE_START
    + 1
)

QHTSF_CHROM = "chr04"
QHTSF_START = 17_867_867
QHTSF_END = 18_681_771


# ================================================================
# Arguments
# ================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Direct 34,649-bp genome-wide matched-window "
            "Hudson FST analysis"
        )
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Input biallelic SNP VCF/VCF.GZ",
    )

    parser.add_argument(
        "--tol",
        required=True,
        help="Tolerant-group sample list",
    )

    parser.add_argument(
        "--sus",
        required=True,
        help="Susceptible-group sample list",
    )

    parser.add_argument(
        "--fai",
        required=True,
        help="Reference FASTA index",
    )

    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix",
    )

    parser.add_argument(
        "--snp-tolerance",
        type=float,
        default=0.20,
        help=(
            "Relative tolerance for informative SNP-count "
            "matching (default: 0.20)"
        ),
    )

    parser.add_argument(
        "--min-informative-snps",
        type=int,
        default=1,
        help=(
            "Minimum informative SNPs required per control "
            "region before matching (default: 1)"
        ),
    )

    return parser.parse_args()


# ================================================================
# Read sample list
# ================================================================

def read_sample_list(path):
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(
            f"Sample-list file not found: {path}"
        )

    samples = []

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            sample = raw.strip()

            if sample and not sample.startswith("#"):
                samples.append(sample)

    if not samples:
        raise ValueError(
            f"Empty sample list: {path}"
        )

    if len(samples) != len(set(samples)):
        duplicates = sorted(
            {
                sample
                for sample in samples
                if samples.count(sample) > 1
            }
        )

        raise ValueError(
            f"Duplicate samples in {path}: {duplicates}"
        )

    return samples


# ================================================================
# Read FAI
# ================================================================

def read_fai(path):
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(
            f"FAI not found: {path}"
        )

    chromosomes = {}
    chromosome_order = []

    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(
            handle,
            start=1,
        ):
            fields = raw.rstrip("\n").split("\t")

            if len(fields) < 2:
                raise ValueError(
                    f"Invalid FAI line {line_number}: "
                    f"{raw.rstrip()}"
                )

            chrom = fields[0].strip()

            try:
                length = int(fields[1])
            except ValueError as error:
                raise ValueError(
                    f"Invalid FAI chromosome length at "
                    f"line {line_number}"
                ) from error

            if chrom in chromosomes:
                raise ValueError(
                    f"Duplicate FAI chromosome: {chrom}"
                )

            if length <= 0:
                raise ValueError(
                    f"Invalid chromosome length: "
                    f"{chrom}={length}"
                )

            chromosomes[chrom] = length
            chromosome_order.append(chrom)

    if CANDIDATE_CHROM not in chromosomes:
        raise ValueError(
            f"{CANDIDATE_CHROM} is absent from the FAI"
        )

    if chromosomes[CANDIDATE_CHROM] < CANDIDATE_END:
        raise ValueError(
            "Candidate interval extends beyond the "
            "chromosome length"
        )

    return chromosomes, chromosome_order


# ================================================================
# Read VCF
# ================================================================

def read_vcf(path):
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(
            f"VCF not found: {path}"
        )

    print("[INFO] Reading VCF...", flush=True)

    callset = allel.read_vcf(
        str(path),
        fields=[
            "samples",
            "variants/CHROM",
            "variants/POS",
            "calldata/GT",
        ],
    )

    if callset is None:
        raise ValueError(
            "No data were read from the VCF"
        )

    required = {
        "samples",
        "variants/CHROM",
        "variants/POS",
        "calldata/GT",
    }

    missing = required - set(callset)

    if missing:
        raise ValueError(
            f"VCF fields not found: {sorted(missing)}"
        )

    samples = np.asarray(
        callset["samples"]
    ).astype(str)

    chrom = np.asarray(
        callset["variants/CHROM"]
    ).astype(str)

    pos = np.asarray(
        callset["variants/POS"],
        dtype=np.int64,
    )

    gt = allel.GenotypeArray(
        callset["calldata/GT"]
    )

    return samples, chrom, pos, gt


# ================================================================
# Validation
# ================================================================

def validate_samples(
    vcf_samples,
    tolerant,
    susceptible,
):
    vcf_set = set(vcf_samples)
    tolerant_set = set(tolerant)
    susceptible_set = set(susceptible)

    overlap = tolerant_set & susceptible_set

    if overlap:
        raise ValueError(
            "Samples occur in both groups: "
            + ", ".join(sorted(overlap))
        )

    missing_tolerant = tolerant_set - vcf_set
    missing_susceptible = susceptible_set - vcf_set

    if missing_tolerant:
        raise ValueError(
            "Tolerant samples absent from VCF: "
            + ", ".join(sorted(missing_tolerant))
        )

    if missing_susceptible:
        raise ValueError(
            "Susceptible samples absent from VCF: "
            + ", ".join(sorted(missing_susceptible))
        )

    sample_index = {
        sample: index
        for index, sample in enumerate(vcf_samples)
    }

    tolerant_indices = np.asarray(
        [sample_index[sample] for sample in tolerant],
        dtype=int,
    )

    susceptible_indices = np.asarray(
        [sample_index[sample] for sample in susceptible],
        dtype=int,
    )

    print("[INPUT VALIDATION 1]")
    print(
        f"  VCF samples:             "
        f"{len(vcf_samples)}"
    )
    print(
        f"  Tolerant samples:        "
        f"{len(tolerant_indices)}"
    )
    print(
        f"  Susceptible samples:     "
        f"{len(susceptible_indices)}"
    )
    print(
        "  Tolerant:                "
        + ", ".join(tolerant)
    )
    print(
        "  Susceptible:             "
        + ", ".join(susceptible)
    )
    print("  Sample validation:       PASS")

    return tolerant_indices, susceptible_indices


def validate_variants(
    chrom,
    pos,
    gt,
    chromosomes,
):
    if len(chrom) != len(pos):
        raise ValueError(
            "CHROM and POS lengths differ"
        )

    if gt.shape[0] != len(pos):
        raise ValueError(
            "Genotype and position lengths differ"
        )

    if np.any(pos < 1):
        raise ValueError(
            "VCF contains POS < 1"
        )

    observed_chromosomes = set(chrom)
    fai_chromosomes = set(chromosomes)

    vcf_only = observed_chromosomes - fai_chromosomes

    if vcf_only:
        raise ValueError(
            "VCF chromosomes absent from FAI: "
            + ", ".join(sorted(vcf_only))
        )

    # -1 is missing; alleles >1 indicate multiallelic genotypes.
    called_genotypes = gt.values[
        gt.values >= 0
    ]

    if (
        called_genotypes.size > 0
        and np.max(called_genotypes) > 1
    ):
        raise ValueError(
            "Multiallelic genotypes were detected. "
            "Use a normalized biallelic SNP-only VCF."
        )

    print("[INPUT VALIDATION 2]")
    print(
        f"  VCF variants:            "
        f"{len(pos):,}"
    )
    print(
        f"  VCF chromosomes:         "
        f"{len(observed_chromosomes)}"
    )
    print(
        f"  Genotype array shape:    "
        f"{gt.shape}"
    )
    print("  Coordinate validation:   PASS")
    print("  Biallelic validation:    PASS")


def validate_sorting(chrom, pos):
    frame = pd.DataFrame(
        {
            "CHROM": chrom,
            "POS": pos,
        }
    )

    duplicate_count = int(
        frame.duplicated(
            subset=["CHROM", "POS"]
        ).sum()
    )

    unsorted_chromosomes = []

    for chromosome, group in frame.groupby(
        "CHROM",
        sort=False,
    ):
        positions = group["POS"].to_numpy()

        if np.any(np.diff(positions) < 0):
            unsorted_chromosomes.append(
                chromosome
            )

    if unsorted_chromosomes:
        raise ValueError(
            "VCF positions are not sorted in: "
            + ", ".join(unsorted_chromosomes)
        )

    print("[INPUT VALIDATION 3]")
    print(
        f"  Duplicate CHROM:POS:     "
        f"{duplicate_count:,}"
    )
    print("  Within-chromosome order: PASS")

    if duplicate_count > 0:
        print(
            "  [WARNING] Duplicate positions are present. "
            "This is acceptable only when they represent "
            "distinct normalized alleles."
        )


# ================================================================
# Hudson FST components
# ================================================================

def calculate_hudson_components(
    gt,
    tolerant_indices,
    susceptible_indices,
):
    print(
        "[INFO] Calculating per-variant Hudson FST "
        "components...",
        flush=True,
    )

    tolerant_gt = gt[:, tolerant_indices]
    susceptible_gt = gt[:, susceptible_indices]

    tolerant_ac = tolerant_gt.count_alleles(
        max_allele=1
    )

    susceptible_ac = susceptible_gt.count_alleles(
        max_allele=1
    )

    numerator, denominator = allel.hudson_fst(
        tolerant_ac,
        susceptible_ac,
    )

    numerator = np.asarray(
        numerator,
        dtype=float,
    )

    denominator = np.asarray(
        denominator,
        dtype=float,
    )

    informative = (
        np.isfinite(numerator)
        & np.isfinite(denominator)
        & (denominator > 0)
    )

    print(
        f"[INFO] Informative variants: "
        f"{np.count_nonzero(informative):,} / "
        f"{len(informative):,}"
    )

    return numerator, denominator, informative


# ================================================================
# Candidate FST
# ================================================================

def calculate_candidate(
    chrom,
    pos,
    numerator,
    denominator,
    informative,
):
    candidate_mask = (
        (chrom == CANDIDATE_CHROM)
        & (pos >= CANDIDATE_START)
        & (pos <= CANDIDATE_END)
        & informative
    )

    n_variants = int(
        np.count_nonzero(candidate_mask)
    )

    if n_variants == 0:
        raise ValueError(
            "No informative variants were found in "
            "the candidate interval"
        )

    sum_numerator = float(
        np.sum(numerator[candidate_mask])
    )

    sum_denominator = float(
        np.sum(denominator[candidate_mask])
    )

    if (
        not np.isfinite(sum_denominator)
        or sum_denominator <= 0
    ):
        raise ValueError(
            "Candidate denominator is not positive"
        )

    candidate_fst = (
        sum_numerator
        / sum_denominator
    )

    result = {
        "CHROM": CANDIDATE_CHROM,
        "START": CANDIDATE_START,
        "END": CANDIDATE_END,
        "LENGTH": REGION_LENGTH,
        "N_INFORMATIVE_VARIANTS": n_variants,
        "SUM_NUMERATOR": sum_numerator,
        "SUM_DENOMINATOR": sum_denominator,
        "HUDSON_FST": candidate_fst,
    }

    print("[CANDIDATE RESULT]")
    print(
        f"  Interval:                "
        f"{CANDIDATE_CHROM}:"
        f"{CANDIDATE_START:,}-"
        f"{CANDIDATE_END:,}"
    )
    print(
        f"  Length:                  "
        f"{REGION_LENGTH:,} bp"
    )
    print(
        f"  Informative variants:    "
        f"{n_variants:,}"
    )
    print(
        f"  Sum numerator:           "
        f"{sum_numerator:.10g}"
    )
    print(
        f"  Sum denominator:         "
        f"{sum_denominator:.10g}"
    )
    print(
        f"  Regional Hudson FST:     "
        f"{candidate_fst:.10g}"
    )

    return result


# ================================================================
# Genome-wide controls
# ================================================================

def overlaps_qhtsf(chromosome, start, end):
    if chromosome != QHTSF_CHROM:
        return False

    return (
        start <= QHTSF_END
        and end >= QHTSF_START
    )


def calculate_control_regions(
    chrom,
    pos,
    numerator,
    denominator,
    informative,
    chromosomes,
    chromosome_order,
    min_informative_snps,
):
    records = []

    print(
        "[INFO] Calculating direct 34,649-bp "
        "genome-wide control FST...",
        flush=True,
    )

    for chromosome in chromosome_order:
        chromosome_length = chromosomes[chromosome]

        if chromosome_length < REGION_LENGTH:
            continue

        chromosome_mask = (
            (chrom == chromosome)
            & informative
        )

        chromosome_pos = pos[
            chromosome_mask
        ]

        chromosome_num = numerator[
            chromosome_mask
        ]

        chromosome_den = denominator[
            chromosome_mask
        ]

        # Only complete non-overlapping regions are retained.
        n_complete_regions = (
            chromosome_length
            // REGION_LENGTH
        )

        if n_complete_regions == 0:
            continue

        if len(chromosome_pos) == 0:
            continue

        # 0-based region index from 1-based VCF position.
        region_index = (
            (chromosome_pos - 1)
            // REGION_LENGTH
        )

        within_complete_region = (
            region_index < n_complete_regions
        )

        region_index = region_index[
            within_complete_region
        ]

        chromosome_num = chromosome_num[
            within_complete_region
        ]

        chromosome_den = chromosome_den[
            within_complete_region
        ]

        counts = np.bincount(
            region_index,
            minlength=n_complete_regions,
        )

        numerator_sums = np.bincount(
            region_index,
            weights=chromosome_num,
            minlength=n_complete_regions,
        )

        denominator_sums = np.bincount(
            region_index,
            weights=chromosome_den,
            minlength=n_complete_regions,
        )

        for index in range(n_complete_regions):
            region_start = (
                index * REGION_LENGTH
                + 1
            )

            region_end = (
                region_start
                + REGION_LENGTH
                - 1
            )

            if overlaps_qhtsf(
                chromosome,
                region_start,
                region_end,
            ):
                continue

            n_variants = int(
                counts[index]
            )

            if n_variants < min_informative_snps:
                continue

            sum_denominator = float(
                denominator_sums[index]
            )

            if (
                not np.isfinite(sum_denominator)
                or sum_denominator <= 0
            ):
                continue

            sum_numerator = float(
                numerator_sums[index]
            )

            region_fst = (
                sum_numerator
                / sum_denominator
            )

            if not np.isfinite(region_fst):
                continue

            records.append(
                {
                    "CHROM": chromosome,
                    "START": region_start,
                    "END": region_end,
                    "LENGTH": REGION_LENGTH,
                    "N_INFORMATIVE_VARIANTS": (
                        n_variants
                    ),
                    "SUM_NUMERATOR": (
                        sum_numerator
                    ),
                    "SUM_DENOMINATOR": (
                        sum_denominator
                    ),
                    "HUDSON_FST": (
                        region_fst
                    ),
                }
            )

    if not records:
        raise ValueError(
            "No valid genome-wide control regions "
            "were generated"
        )

    controls = pd.DataFrame(records)

    print(
        f"[INFO] Valid control regions: "
        f"{len(controls):,}"
    )

    return controls


# ================================================================
# Match controls and calculate empirical P
# ================================================================

def match_and_test(
    controls,
    candidate,
    snp_tolerance,
):
    candidate_n = candidate[
        "N_INFORMATIVE_VARIANTS"
    ]

    lower = candidate_n * (
        1 - snp_tolerance
    )

    upper = candidate_n * (
        1 + snp_tolerance
    )

    matched = controls[
        controls[
            "N_INFORMATIVE_VARIANTS"
        ].between(
            lower,
            upper,
            inclusive="both",
        )
    ].copy()

    if matched.empty:
        raise ValueError(
            "No SNP-count-matched controls remained. "
            "Increase --snp-tolerance."
        )

    candidate_fst = candidate[
        "HUDSON_FST"
    ]

    n_matched = len(matched)

    n_equal_or_higher = int(
        np.count_nonzero(
            matched["HUDSON_FST"].to_numpy()
            >= candidate_fst
        )
    )

    n_strictly_higher = int(
        np.count_nonzero(
            matched["HUDSON_FST"].to_numpy()
            > candidate_fst
        )
    )

    empirical_p = (
        n_equal_or_higher + 1
    ) / (
        n_matched + 1
    )

    empirical_percentile = float(
        np.mean(
            matched["HUDSON_FST"].to_numpy()
            < candidate_fst
        )
        * 100
    )

    matched["MATCH_LOWER_N_VARIANTS"] = (
        lower
    )

    matched["MATCH_UPPER_N_VARIANTS"] = (
        upper
    )

    matched["CANDIDATE_N_VARIANTS"] = (
        candidate_n
    )

    matched["CANDIDATE_HUDSON_FST"] = (
        candidate_fst
    )

    matched["EQUAL_OR_HIGHER_THAN_CANDIDATE"] = (
        matched["HUDSON_FST"]
        >= candidate_fst
    )

    matched["EMPIRICAL_P"] = (
        empirical_p
    )

    statistics = {
        "SNP_MATCH_LOWER": lower,
        "SNP_MATCH_UPPER": upper,
        "N_MATCHED_CONTROLS": n_matched,
        "N_EQUAL_OR_HIGHER": n_equal_or_higher,
        "N_STRICTLY_HIGHER": n_strictly_higher,
        "EMPIRICAL_P": empirical_p,
        "EMPIRICAL_PERCENTILE": (
            empirical_percentile
        ),
    }

    return matched, statistics


# ================================================================
# Main
# ================================================================

def main():
    args = parse_args()

    if not 0 <= args.snp_tolerance < 1:
        raise ValueError(
            "--snp-tolerance must be >= 0 and < 1"
        )

    if args.min_informative_snps < 1:
        raise ValueError(
            "--min-informative-snps must be >= 1"
        )

    tolerant = read_sample_list(
        args.tol
    )

    susceptible = read_sample_list(
        args.sus
    )

    chromosomes, chromosome_order = read_fai(
        args.fai
    )

    samples, chrom, pos, gt = read_vcf(
        args.vcf
    )

    tolerant_indices, susceptible_indices = (
        validate_samples(
            samples,
            tolerant,
            susceptible,
        )
    )

    validate_variants(
        chrom,
        pos,
        gt,
        chromosomes,
    )

    validate_sorting(
        chrom,
        pos,
    )

    numerator, denominator, informative = (
        calculate_hudson_components(
            gt,
            tolerant_indices,
            susceptible_indices,
        )
    )

    candidate = calculate_candidate(
        chrom,
        pos,
        numerator,
        denominator,
        informative,
    )

    controls = calculate_control_regions(
        chrom,
        pos,
        numerator,
        denominator,
        informative,
        chromosomes,
        chromosome_order,
        args.min_informative_snps,
    )

    matched, statistics = match_and_test(
        controls,
        candidate,
        args.snp_tolerance,
    )

    prefix = Path(
        args.out_prefix
    )

    all_controls_path = Path(
        f"{prefix}.all_controls.tsv"
    )

    matched_path = Path(
        f"{prefix}.matched_controls.tsv"
    )

    candidate_path = Path(
        f"{prefix}.candidate.tsv"
    )

    summary_path = Path(
        f"{prefix}.summary.tsv"
    )

    controls = controls.sort_values(
        ["CHROM", "START"]
    )

    matched = matched.sort_values(
        ["HUDSON_FST", "CHROM", "START"],
        ascending=[False, True, True],
    )

    controls.to_csv(
        all_controls_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    matched.to_csv(
        matched_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    pd.DataFrame(
        [candidate]
    ).to_csv(
        candidate_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary = {
        **candidate,
        "SNP_TOLERANCE": (
            args.snp_tolerance
        ),
        "TOTAL_CONTROL_REGIONS": (
            len(controls)
        ),
        **statistics,
    }

    pd.DataFrame(
        [summary]
    ).to_csv(
        summary_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print("[FINAL RESULTS]")
    print(
        f"  Candidate FST:           "
        f"{candidate['HUDSON_FST']:.10g}"
    )
    print(
        f"  Candidate SNPs:          "
        f"{candidate['N_INFORMATIVE_VARIANTS']:,}"
    )
    print(
        f"  SNP matching range:      "
        f"{statistics['SNP_MATCH_LOWER']:.2f}-"
        f"{statistics['SNP_MATCH_UPPER']:.2f}"
    )
    print(
        f"  Total controls:          "
        f"{len(controls):,}"
    )
    print(
        f"  Matched controls:        "
        f"{statistics['N_MATCHED_CONTROLS']:,}"
    )
    print(
        f"  Controls >= candidate:   "
        f"{statistics['N_EQUAL_OR_HIGHER']:,}"
    )
    print(
        f"  Controls > candidate:    "
        f"{statistics['N_STRICTLY_HIGHER']:,}"
    )
    print(
        f"  Empirical P:             "
        f"{statistics['EMPIRICAL_P']:.10g}"
    )
    print(
        f"  Empirical percentile:    "
        f"{statistics['EMPIRICAL_PERCENTILE']:.4f}"
    )
    print(
        f"  Candidate output:        "
        f"{candidate_path}"
    )
    print(
        f"  All controls output:     "
        f"{all_controls_path}"
    )
    print(
        f"  Matched controls output: "
        f"{matched_path}"
    )
    print(
        f"  Summary output:          "
        f"{summary_path}"
    )

    print(
        "\n[INTERPRETATION]\n"
        "  The empirical P value measures how often a "
        "genome-wide 34,649-bp\n"
        "  region with a similar number of informative "
        "SNPs has Hudson FST equal\n"
        "  to or greater than the candidate interval. "
        "It is not a phenotype-\n"
        "  association P value and does not correct for "
        "population structure."
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(
            f"[ERROR] {error}",
            file=sys.stderr,
        )
        sys.exit(1)
