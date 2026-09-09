# Structural and sequence-level dissection of qHTSF4.1 reveals a chromosomal inversion and a candidate interval in rice

## Overview

 This repository contains scripts and derived datasets supporting structural and population-genomic analyses of the rice flowering-stage heat-tolerance locus qHTSF4.1. The deposited materials cover genome-wide Hudson's FST analysis, a sensitivity analysis excluding Giza178, sequence-identity-based introgression mapping, visualization of introgressed regions and previously reported reproductive-stage heat-tolerance QTLs, and long-read alignment analysis at the corresponding inversion boundaries.

## Repository contents

| File | Description |
|---|---|
| `genomewide_hudson_fst.py` | Calculates Hudson's *F*ST from callable biallelic SNPs and summarizes differentiation in non-overlapping 10-kb windows. |
| `hudson_fst_5samples_10kb_nonoverlapping.tsv` | Genome-wide 10-kb window results for the comparison of N22, Giza178, and IR64-HTSF4.1 with Azucena and Nipponbare. |
| `hudson_fst_4samples_noGiza178_10kb_nonoverlapping.tsv` | Sensitivity-analysis results for the comparison of N22 and IR64-HTSF4.1 with Azucena and Nipponbare after excluding Giza178. |
| `qHTSF4.1_5samples_callable_biallelic.vcf` | Callable biallelic variants within the *qHTSF4.1* region for N22, Giza178, IR64-HTSF4.1, Azucena, and Nipponbare. |
| `detect_introgression_regions_nucmer.R` | Processes sequence-identity information obtained from NUCmer and `show-coords` alignments and assigns IR64-HTSF4.1 genomic windows to the more similar parental reference. |
| `plot_heatQTL_introgression_regions.py` | Visualizes N22-assigned genomic segments together with previously reported N22-derived reproductive-stage heat-tolerance QTLs and linked markers after coordinate conversion to the IR64-HTSF4.1 assembly. |
| `analyze_breakpoint_read_support.py` | Counts continuous, clipped, and supplementary ONT read alignments at the corresponding inversion boundaries in the N22 and IR64sj assemblies. |

## Expected outputs and key validation results

The analyses generated the following principal reproducibility checkpoints:

- Whole-genome alignment identified an approximately 810.3-kb inversion corresponding to Chr04:17,207,988-18,018,255 in the IR64-HTSF4.1 assembly.
- The LD-defined candidate interval was Chr04:18,535,070-18,569,718 (34,649 bp) in the Nipponbare reference assembly (GCA_001433935.1).
- Regional Hudson's *F*ST across the candidate interval was 0.945. Six of 192 genomic control regions matched for interval length and SNP count had regional *F*ST values greater than or equal to 0.945.
- Elevated differentiation within the candidate interval was retained in the sensitivity analysis excluding Giza178.
- IR64-HTSF4.1 ONT reads produced 75 and 78 continuous alignments spanning the left and right corresponding boundaries, respectively, in the N22 assembly. No continuous alignment spanned either corresponding boundary in the IR64sj assembly, where boundary-associated clipping and supplementary alignments were observed.

These results provide validation targets for reproducing the deposited analyses. They do not establish the inversion, candidate interval, or individual variants as causal determinants of flowering-stage heat tolerance.

## Plant materials and reference assemblies

IR64-HTSF4.1 is an IR64-derived backcross line developed using Nagina 22 (N22) as the donor of the *qHTSF4.1* region. The line analyzed in the genome-assembly work was derived from BC2F7 plants. The available material identifiers are IR64 accession IRGC 117268 (GID 4537782) and IR64-HTSF4.1 GID 4537748.

| Genome or material | Role in the analyses | Accession or identifier |
|---|---|---|
| IR64sj | IR64 reference used for scaffolding, introgression mapping, structural comparison, and long-read alignment | GCA_046630055.1 |
| Nagina 22 (N22) | Donor reference used for introgression mapping and long-read alignment | GCA_001952365.2 |
| Nipponbare | Reference used for five-accession variant calling and candidate-interval coordinates | GCA_001433935.1 |
| Azucena | Heat-susceptible accession included in the five-accession comparison | GCA_009830595.1 |
| IR64-HTSF4.1 | N22-derived *qHTSF4.1* introgression line in the IR64 genetic background | GID 4537748 |

## Five-accession Hudson's *F*ST analysis

Accessions were grouped according to flowering-stage heat-tolerance classifications reported in previous studies:

| Analysis group | Accessions |
|---|---|
| Previously reported as heat tolerant | N22, Giza178, and IR64-HTSF4.1 |
| Previously reported as heat susceptible | Azucena and Nipponbare |

The input VCF was filtered to retain callable biallelic SNPs with non-missing diploid genotypes in every accession included in the comparison. Hudson's *F*ST was calculated using scikit-allel v1.3.13 and summarized across non-overlapping 10-kb genomic windows. The 90th percentile of valid genome-wide windows was used as the empirical background threshold.

The Giza178-exclusion sensitivity analysis compared N22 and IR64-HTSF4.1 with Azucena and Nipponbare using the same non-overlapping 10-kb window framework.

The heat-tolerance classifications were obtained from previous studies rather than contemporaneous phenotyping under a single experimental protocol. The *F*ST results should therefore be interpreted as exploratory differentiation signals rather than phenotype-aware association results.

## Regional VCF and candidate interval

`qHTSF4.1_5samples_callable_biallelic.vcf` contains callable biallelic variants for all five accessions within the broader *qHTSF4.1* region. Variants were called relative to the Nipponbare reference assembly (GCA_001433935.1).

The LD-defined candidate interval evaluated in the manuscript is:

```text
Chr04:18,535,070-18,569,718
Length: 34,649 bp
Reference: Nipponbare, GCA_001433935.1
```

The exact boundaries of the broader *qHTSF4.1* region represented in the VCF should be obtained from the VCF records and documented here before archival release.

Genome-wide differentiation was evaluated in non-overlapping 10-kb windows. Two of the four windows overlapping the 34,649-bp candidate interval exceeded the genome-wide 90th percentile:

| Window | Hudson's *F*ST | Informative SNPs |
|---|---:|---:|
| Chr04:18,531,610-18,541,609 | 0.908333 | 18 |
| Chr04:18,541,610-18,551,609 | 0.922581 | 12 |

The regional Hudson's *F*ST across the complete 34,649-bp interval was 0.945. For the matched-region analysis, non-overlapping genomic intervals of identical length and with SNP counts within +/-20% of the candidate-interval count were used as controls after excluding intervals overlapping the broader *qHTSF4.1* region. Six of 192 matched control regions had regional *F*ST values greater than or equal to 0.945.

This five-accession VCF was used for regional *F*ST and candidate-variant analyses. It was not the population-level LD input for the analyses of the complete 3,000 Rice Genomes panel and its subpopulations.

## Introgression mapping

The chromosome-level IR64-HTSF4.1 assembly was divided into non-overlapping 10-kb windows, and the corresponding sequences were extracted using SAMtools. Each window was aligned independently to the IR64sj and N22 reference assemblies using NUCmer from MUMmer4. Alignment coordinates and percentage sequence identities were obtained using `show-coords`.

Each IR64-HTSF4.1 window was assigned to the parental reference with the higher sequence identity. Adjacent N22-assigned windows were interpreted collectively as N22-associated introgressed regions. Because ancestry assignment was based on fixed 10-kb windows, the inferred boundaries and chromosome-wide introgressed proportions are approximate.

