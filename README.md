# Structural and sequence-level dissection of qHTSF4.1 reveals a chromosomal inversion and a candidate interval in rice

## Overview

 This repository contains scripts and derived datasets supporting structural and population-genomic analyses of the rice flowering-stage heat-tolerance locus qHTSF4.1. The deposited materials cover genome-wide Hudson's FST analysis, a sensitivity analysis excluding Giza178, sequence-identity-based introgression mapping, visualization of introgressed regions and previously reported reproductive-stage heat-tolerance QTLs, and long-read alignment analysis at the corresponding inversion boundaries.

## Repository contents

| File | Purpose |
|---|---|
| `genomewide_hudson_fst.py` | Calculates Hudson's *F*ST from callable biallelic SNPs and summarizes differentiation in non-overlapping 10-kb windows. |
| `hudson_fst_5samples_10kb_nonoverlapping.tsv` | Contains genome-wide 10-kb window results for the five-accession comparison. |
| `hudson_fst_4samples_noGiza178_10kb_nonoverlapping.tsv` | Contains 10-kb window results from the sensitivity analysis excluding Giza178. |
| `qHTSF4.1_5samples_callable_biallelic.vcf` | Contains callable biallelic variants within *qHTSF4.1* for the five accessions used in the regional analysis. |
| `detect_introgression_regions_nucmer.R` | Uses sequence-identity information derived from NUCmer and `show-coords` alignments to assign IR64-HTSF4.1 genomic windows to the more similar parental reference. |
| `plot_heatQTL_introgression_regions.py` | Plots N22-assigned genomic segments and previously reported N22-derived reproductive-stage heat-tolerance QTLs or linked markers in the IR64-HTSF4.1 coordinate system. |
| `analyze_breakpoint_read_support.py` | Classifies ONT read alignments as continuous, clipped, or clipped with supplementary-alignment information at configured inversion boundaries. |

## Plant materials and reference assemblies

IR64-HTSF4.1 is an IR64-derived backcross line developed using Nagina 22 (N22) as the donor of the *qHTSF4.1* region. The genome assembly and annotation were generated from BC2F7 plants.

| Genome or material | Role | Accession or identifier |
|---|---|---|
| IR64sj | IR64 reference used for scaffolding, introgression mapping, structural comparison, and long-read alignment | GCA_046630055.1 |
| Nagina 22 (N22) | Donor reference used for introgression mapping and long-read alignment | GCA_001952365.2 |
| Nipponbare | Reference used for five-accession variant calling and candidate-region coordinates | GCA_001433935.1 |
| Azucena | Accession included in the five-accession comparison | GCA_009830595.1 |
| IR64 | Recurrent parent | IRGC 117268; GID 4537782 |
| IR64-HTSF4.1 | N22-derived *qHTSF4.1* introgression line in the IR64 genetic background | GID 4537748 |

## Analysis 1: Hudson's *F*ST

### Purpose

This analysis evaluates sequence differentiation between accessions previously classified as heat tolerant and heat susceptible. Genome-wide estimates provide an empirical background for interpreting differentiation within *qHTSF4.1*. A second analysis excluding Giza178 evaluates whether the regional profile depends on its inclusion in the tolerant group.

### Sample groups

| Analysis | Previously reported as heat tolerant | Previously reported as heat susceptible |
|---|---|---|
| Five-accession analysis | N22, Giza178, IR64-HTSF4.1 | Azucena, Nipponbare |
| Sensitivity analysis | N22, IR64-HTSF4.1 | Azucena, Nipponbare |

The classifications were obtained from previous studies rather than phenotyping conducted concurrently with the genomic analyses.

### Input and filtering

The genome-wide analysis requires a multi-sample VCF containing the five accessions. The VCF was filtered to retain callable biallelic SNPs with non-missing diploid genotypes in all accessions included in each comparison. Hudson's *F*ST was calculated using scikit-allel v1.3.13 and summarized in non-overlapping 10-kb windows.

`qHTSF4.1_5samples_callable_biallelic.vcf` is the regional five-accession VCF used for analysis of variants within *qHTSF4.1*. It is not a substitute for the genome-wide VCF required to regenerate the genome-wide window results.

The LD-defined candidate interval examined within the regional VCF is Chr04:18,535,070-18,569,718 in the Nipponbare reference assembly (GCA_001433935.1).

## Analysis 2: Introgression-region assignment

### Purpose

This analysis identifies genomic windows in IR64-HTSF4.1 that are more similar to the N22 donor reference than to the IR64 recurrent-parent reference. The resulting assignments are used to describe the genome-wide distribution of donor-associated segments.

### Inputs

- IR64-HTSF4.1 chromosome-level genome assembly
- IR64sj reference assembly (GCA_046630055.1)
- N22 reference assembly (GCA_001952365.2)
- NUCmer alignments of IR64-HTSF4.1 windows to each parental reference
- Percentage sequence identities and alignment coordinates obtained with `show-coords`

The IR64-HTSF4.1 assembly was divided into non-overlapping 10-kb windows, and each window was aligned independently to IR64sj and N22 using NUCmer from MUMmer4. `detect_introgression_regions_nucmer.R` assigns each window to the parental reference with the higher sequence identity. Because the analysis uses fixed 10-kb windows, inferred ancestry boundaries are approximate.

## Analysis 3: Heat-QTL and introgression visualization

### Purpose

This analysis compares N22-assigned genomic segments with previously reported N22-derived reproductive-stage heat-tolerance QTLs and linked markers.

### Inputs

- Introgression-region assignments generated from the window-based sequence-identity analysis
- Published QTL intervals and linked-marker positions collected from Ye et al. (2012), Shanmugavadivel et al. (2017), and Nguyen et al. (2022)
- QTL and marker positions converted to the IR64-HTSF4.1 coordinate system
- IR64-HTSF4.1 chromosome lengths

Published QTL and marker coordinates were converted to the IR64-HTSF4.1 assembly using minimap2 and BLASTN sequence alignments. `plot_heatQTL_introgression_regions.py` jointly visualizes the converted loci and inferred N22-assigned segments.

## Analysis 4: Long-read support at inversion boundaries

### Purpose

This analysis examines the original IR64-HTSF4.1 ONT reads independently of the scaffolded IR64-HTSF4.1 pseudomolecules to evaluate read-alignment patterns at the corresponding inversion boundaries in N22 and IR64sj.

### Inputs

- Original IR64-HTSF4.1 ONT reads
- N22 reference assembly (GCA_001952365.2)
- IR64sj reference assembly (GCA_046630055.1)
- Coordinate-sorted and indexed alignment files generated by aligning the ONT reads independently to each reference
- Configured boundary coordinates shown below

| Reference assembly | Boundary | Coordinate (bp) |
|---|---|---:|
| N22 | Left | 16,889,637 |
| N22 | Right | 17,700,080 |
| IR64sj | Left | 17,273,104 |
| IR64sj | Right | 18,073,176 |

ONT reads were aligned with minimap2 v2.30-r1287 using the `map-ont` preset. `analyze_breakpoint_read_support.py` counts single alignments that traverse a configured boundary, reads clipped near a boundary, and clipped reads carrying an `SA` tag.
