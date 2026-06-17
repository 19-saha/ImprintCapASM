# ImprintCapASM <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/19-saha/ImprintCapASM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/19-saha/ImprintCapASM/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/ImprintCapASM)](https://CRAN.R-project.org/package=ImprintCapASM)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview

**ImprintCapASM** is an R package for **SNP-phased allele-specific methylation (ASM) analysis** across the 41 known human imprinted differentially methylated regions (DMRs). It is designed for clinical diagnostic workflows that profile imprint disorder cases — including Beckwith-Wiedemann syndrome (BWS), Silver-Russell syndrome (SRS), Prader-Willi syndrome (PWS), Angelman syndrome (AS), and related conditions — from bisulfite sequencing data produced by targeted capture panels.

The package provides three core functions that form a sequential pipeline:

1. **`prepare_cpg_snp_input()`** — Links CpG methylation values to nearby heterozygous SNPs; produces a per-sample Excel table and a BED file
2. **`extract_bam_regions()`** — Extracts and sorts a BAM subset covering the SNP windows for each sample
3. **`ASM()`** — Reads the extracted BAM, assigns each read to a parental allele, and computes allele-specific methylation statistics; returns three output tables and a line-plot PDF

For processing **multiple samples together** — the standard diagnostic use case — use `run_pipeline()`, which runs the full three-step pipeline for all control samples as a batch, and separately for all patient samples as a batch. Controls and patients are always run independently using their respective `filter_cpgs` reference files.

---

## Background

Genomic imprinting is an epigenetic phenomenon whereby a subset of genes are expressed in a parent-of-origin dependent manner, regulated by differentially methylated regions (DMRs). Loss or gain of methylation at these DMRs underlies a class of rare congenital disorders collectively known as imprinting disorders. Accurate diagnosis requires quantifying the methylation of each parental allele **separately** — a task that standard bisulfite sequencing alone cannot achieve without phasing methylation data to nearby heterozygous SNPs.

ImprintCapASM implements a SNP-phasing strategy: heterozygous SNPs detected in bisulfite sequencing reads are used to assign each read to a parental allele (REF or ALT), and CpG methylation values on each allele are computed and compared. Deviation from the expected allele-specific methylation pattern at a given DMR indicates a potential imprinting disorder.

---

## Installation

### From CRAN (stable release)

```r
install.packages("ImprintCapASM")
```

### From GitHub (development version)

```r
# install.packages("remotes")
remotes::install_github("19-saha/ImprintCapASM")
```

### Bioconductor dependencies

ImprintCapASM depends on several Bioconductor packages. Install them first if not already present:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "BiocParallel",
    "Rsamtools",
    "GenomicRanges",
    "IRanges",
    "S4Vectors",
    "SummarizedExperiment",
    "VariantAnnotation"
))
```

---

## The Two `filter_cpgs` Reference Files

A key concept in ImprintCapASM is that **controls and patients each have their own `filter_cpgs` reference file**. These are not interchangeable:

| File | Used with | Contains | Purpose |
|------|-----------|----------|---------|
| `inst/extdata/filter_cpgs_ctrl.xlsx` | `sample_type = "control"` | `Control_1`, `Control_2`, ... columns | Computes mean/SD methylation and CpG variance categories from the control cohort |
| `inst/extdata/filter_cpgs_pat.xlsx` | `sample_type = "patient"` | `Patient_1`, `Patient_2`, ... columns | Computes mean/SD methylation and CpG variance categories from the patient cohort |

Both files share the same structure (`chr`, `5_location`, `3_location`, `DMR`, then sample columns). The `ASM()` function auto-detects sample columns by matching the pattern `^Control_` or `^Patient_` in the column names. Passing the wrong file to the wrong `sample_type` will produce incorrect variance categories and misleading plots.

Both files are used identically for CpG window definition in `prepare_cpg_snp_input()` — what differs is the cohort-specific methylation statistics computed during `ASM()`.

---

## Pipeline Overview

```
Bisulfite sequencing run (targeted imprint capture panel)
        │
        ├── bssnper SNP calling    →  sample.SNPs.out      (VCFv4.3, plain text)
        ├── bssnper CG methylation →  sample.CGmeth.txt    (9-column TSV)
        └── Picard MarkDuplicates  →  sample_markdup.bam + .bai
                │
                ▼  [per sample, run separately for controls and patients]
        ┌──────────────────────────────────────────┐
        │   prepare_cpg_snp_input()                │
        │   Input:  sample.SNPs.out                │  Filters heterozygous SNPs (GT=0/1),
        │           sample.CGmeth.txt              │  overlaps with CpG panel windows,
        │           inst/extdata/filter_cpgs_ctrl.xlsx          │  joins CpG methylation fractions
        │             OR inst/extdata/filter_cpgs_pat.xlsx      │
        │   Output: cpg_snps_CG_{type}_{id}.xlsx   │
        │           cpg_snps_CG_{type}_{id}.bed    │
        └──────────────┬───────────────────────────┘
                       │
        ┌──────────────▼───────────────────────────┐
        │   extract_bam_regions()                  │
        │   Input:  sample_markdup.bam             │  Subsets BAM to SNP windows,
        │           cpg_snps_CG_{type}_{id}.bed    │  sorts and indexes the output
        │   Output: {type}_{id}_wide.bam + .bai    │
        └──────────────┬───────────────────────────┘
                       │
             [all samples of same type combined]
                       │
        ┌──────────────▼───────────────────────────┐
        │   ASM()                                  │
        │   Input:  cpg_snps_CG_{type}_{id}.xlsx   │  Bisulfite-aware allele assignment,
        │           {type}_{id}_wide.bam           │  per-read methylation scoring,
        │           inst/extdata/filter_cpgs_ctrl.xlsx          │  CpG variance classification using
        │             OR inst/extdata/filter_cpgs_pat.xlsx      │  cohort-matched reference
        │   Output: asm_{type}_{id}.xlsx           │
        │           snp_cpg_{type}_{id}.xlsx       │
        │           meth_summary_{type}_{id}.xlsx  │
        │           lineplot_{type}_{id}.pdf       │
        └──────────────────────────────────────────┘
```

---

## Input File Formats

### 1. SNP file — `sample.SNPs.out` (bssnper VCFv4.3)

Produced by [BS-Snper](https://github.com/hellbelly/BS-Snper). Plain-text VCF — **no bgzip or tabix index required**. The function reads the `GT` FORMAT field and retains only heterozygous SNPs (`GT == "0/1"`) with sufficient depth:

```
#CHROM  POS     ID  REF ALT QUAL  FILTER  INFO             FORMAT                                 SAMPLE
chr11   2016400 .   G   A   85    PASS    DP=28;AD=15,13;  GT:DP:AD:ADF:ADR:BSD:BSQ:ALFR          0/1:28:15,13:...
```

### 2. Methylation file — `sample.CGmeth.txt` (bssnper CG output)

Tab-delimited, 9 columns with a `#CHROM` header. Watson and Crick strand methylation and coverage are merged internally by the function:

```
#CHROM  POS       CONTEXT  Watson-METH  Watson-COVERAGE  Watson-QUAL  Crick-METH  Crick-COVERAGE  Crick-QUAL
chr11   2016405   CG       155          169              33           365         494             33
```

### 3. CpG panel reference files — `inst/extdata/filter_cpgs_ctrl.xlsx` and `inst/extdata/filter_cpgs_pat.xlsx`

Two separate reference Excel files — one for controls, one for patients. Both share the same column structure: genomic coordinates and DMR name, followed by per-sample methylation percentages. The `ASM()` function detects sample columns automatically by matching `^Control_` or `^Patient_` column name prefixes:

```
chr    5_location  3_location  DMR        Control_1  Control_2  Control_3  ...
chr11  2016404     2016406     H19/IGF2   82         84         81         ...
```

```
chr    5_location  3_location  DMR        Patient_1  Patient_2  Patient_3  ...
chr11  2016404     2016406     H19/IGF2   45         83         80         ...
```

### 4. BAM file — `sample_markdup.bam` + `.bam.bai`

Duplicate-marked, coordinate-sorted BAM produced by Picard `MarkDuplicates`. The `.bai` index must be present alongside the BAM. If the index is missing, `extract_bam_regions()` creates it automatically via `Rsamtools::indexBam()`.

---

## Recommended Folder Structure

Organise your project with controls and patients in separate folders so that `run_pipeline()` can glob files cleanly:

```
project/
├── controls/
│   ├── snps/
│   │   ├── CTRL_01.SNPs.out
│   │   ├── CTRL_02.SNPs.out
│   │   └── ...
│   ├── meth/
│   │   ├── CTRL_01.CGmeth.txt
│   │   ├── CTRL_02.CGmeth.txt
│   │   └── ...
│   ├── bams/
│   │   ├── CTRL_01_markdup.bam
│   │   ├── CTRL_01_markdup.bam.bai
│   │   └── ...
│   └── output/
│
├── patients/
│   ├── snps/
│   ├── meth/
│   ├── bams/
│   └── output/
│
├── inst/extdata/filter_cpgs_ctrl.xlsx   ← control reference panel
└── inst/extdata/filter_cpgs_pat.xlsx    ← patient reference panel
```

---

## Usage

### Running a single control sample

```r
library(ImprintCapASM)

# Step 1
prepare_cpg_snp_input(
    snp_file     = "controls/snps/CTRL_01.SNPs.out",
    meth_file    = "controls/meth/CTRL_01.CGmeth.txt",
    cpg_ref_file = "inst/extdata/filter_cpgs_ctrl.xlsx",
    sample_type  = "control"
)
# Writes: cpg_snps_CG_control_CTRL_01.xlsx
#         cpg_snps_CG_control_CTRL_01.bed

# Step 2
extract_bam_regions(
    bam_file    = "controls/bams/CTRL_01_markdup.bam",
    bed_file    = "cpg_snps_CG_control_CTRL_01.bed",
    output_dir  = "controls/output/",
    sample_type = "control"
)
# Writes: controls/output/control_CTRL_01_wide.bam + .bai

# Step 3
ASM(
    cpg_snp_file     = "cpg_snps_CG_control_CTRL_01.xlsx",
    sam_file         = "controls/output/control_CTRL_01_wide.bam",
    filter_cpgs_file = "inst/extdata/filter_cpgs_ctrl.xlsx",
    sample_type      = "control"
)
# Writes: asm_control_CTRL_01.xlsx
#         snp_cpg_control_CTRL_01.xlsx
#         meth_summary_control_CTRL_01.xlsx
#         lineplot_control_CTRL_01.pdf
```

### Running a single patient sample

```r
# Step 1
prepare_cpg_snp_input(
    snp_file     = "patients/snps/PAT_01.SNPs.out",
    meth_file    = "patients/meth/PAT_01.CGmeth.txt",
    cpg_ref_file = "inst/extdata/filter_cpgs_pat.xlsx",       
    sample_type  = "patient"
)

# Step 2
extract_bam_regions(
    bam_file    = "patients/bams/PAT_01_markdup.bam",
    bed_file    = "cpg_snps_CG_patient_PAT_01.bed",
    output_dir  = "patients/output/",
    sample_type = "patient"
)

# Step 3
ASM(
    cpg_snp_file     = "cpg_snps_CG_patient_PAT_01.xlsx",
    sam_file         = "patients/output/patient_PAT_01_wide.bam",
    filter_cpgs_file = "inst/extdata/filter_cpgs_pat.xlsx",   # <-- patient reference
    sample_type      = "patient"
)
```

---

### Running a full cohort with `run_pipeline()`

`run_pipeline()` processes all samples in a given folder in batch. Controls and patients are always run as **separate calls** with their respective reference files:

```r
library(ImprintCapASM)

# --- Run all controls ---
run_pipeline(
    snp_dir          = "controls/snps/",
    meth_dir         = "controls/meth/",
    bam_dir          = "controls/bams/",
    filter_cpgs_file = "inst/extdata/filter_cpgs_ctrl.xlsx",
    output_dir       = "controls/output/",
    sample_type      = "control"
)

# --- Run all patients (separate call, separate reference file) ---
run_pipeline(
    snp_dir          = "patients/snps/",
    meth_dir         = "patients/meth/",
    bam_dir          = "patients/bams/",
    filter_cpgs_file = "inst/extdata/filter_cpgs_pat.xlsx",
    output_dir       = "patients/output/",
    sample_type      = "patient"
)
```

`run_pipeline()` automatically matches files across `snp_dir`, `meth_dir`, and `bam_dir` by sample ID, iterates Steps 1 and 2 per sample, then calls `ASM()` on the combined output for that cohort.

---

## Toy Example with Built-in Data

The package ships with minimal example files covering two chr11 DMRs (H19/IGF2 and KCNQ1OT1):

```r
library(ImprintCapASM)

snp_file     <- system.file("extdata", "example_snp.vcf",         package = "ImprintCapASM")
meth_file    <- system.file("extdata", "example_cgmeth.txt",       package = "ImprintCapASM")
cpg_ref_file <- system.file("extdata", "example_filter_cpgs.xlsx", package = "ImprintCapASM")
bam_file     <- system.file("extdata", "example.bam",              package = "ImprintCapASM")

# Step 1
prepare_cpg_snp_input(
    snp_file     = snp_file,
    meth_file    = meth_file,
    cpg_ref_file = cpg_ref_file,
    sample_type  = "control"
)

# Step 2
extract_bam_regions(
    bam_file    = bam_file,
    bed_file    = list.files(tempdir(), pattern = "\\.bed$", full.names = TRUE)[1],
    output_dir  = tempdir(),
    sample_type = "control"
)

# Step 3
ASM(
    cpg_snp_file     = list.files(tempdir(), pattern = "cpg_snps.*\\.xlsx$", full.names = TRUE)[1],
    sam_file         = list.files(tempdir(), pattern = "_wide\\.bam$",        full.names = TRUE)[1],
    filter_cpgs_file = cpg_ref_file,
    sample_type      = "control"
)
```

---

## Output Files and Column Descriptions

`ASM()` writes **three Excel files** and **one PDF** per run.

---

### 1. `asm_{type}_{sample_id}.xlsx` — Read-level allele-methylation table

One row per read–CpG combination. The most granular output.

| Column | Description |
|--------|-------------|
| `sample_id` | Sample identifier (derived from BAM filename) |
| `sample_type` | `"control"` or `"patient"` |
| `id` | Read name |
| `read_sequence` | Raw read sequence |
| `read_start` | Leftmost mapping position of the read |
| `flag` | SAM FLAG value |
| `flag_context` | Human-readable FLAG interpretation |
| `combined_tags` | Concatenated SAM optional tags |
| `strand` | `"Forward"` or `"Reverse"` |
| `chr` | Chromosome |
| `DMR` | Imprinted DMR name (e.g. `H19/IGF2`) |
| `snp_pos` | Genomic position of the phasing SNP |
| `cpg_pos` | Genomic position of the CpG |
| `allele_type` | `"REF"` or `"ALT"` (parental allele assignment) |
| `ref_allele` | Reference base at the SNP |
| `alt_allele` | Alternative base at the SNP |
| `assignment_note` | Bisulfite-aware logic used for allele assignment |
| `n_methylated` | 1 if the CpG is methylated on this read, else 0 |
| `n_unmethylated` | 1 if the CpG is unmethylated on this read, else 0 |
| `meth_frac` | Same as `n_methylated` (numeric; used for summaries) |
| `Padded_Sequence` | Read sequence left-padded for DMR alignment visualisation |
| `mean_methylation` | Cohort mean methylation for this CpG (from `filter_cpgs`) |
| `sd_methylation` | Cohort SD for this CpG (from `filter_cpgs`) |
| `Category` | CpG variance class: `LOWvar`, `MSDvar`, `SDvar`, or `Mvar` |

---

### 2. `snp_cpg_{type}_{sample_id}.xlsx` — Per SNP–CpG pair summary

One row per unique (SNP position, CpG position) combination, with allele-stratified read counts and methylation fractions.

| Column | Description |
|--------|-------------|
| `snp_pos` | Genomic position of the phasing SNP |
| `cpg_pos` | Genomic position of the CpG |
| `sample_id` | Sample identifier |
| `chr` | Chromosome |
| `DMR` | Imprinted DMR name |
| `ref_allele` | Reference base at the SNP |
| `alt_allele` | Alternative base at the SNP |
| `REF_m` | Methylated read count on the REF allele |
| `REF_um` | Unmethylated read count on the REF allele |
| `ALT_m` | Methylated read count on the ALT allele |
| `ALT_um` | Unmethylated read count on the ALT allele |
| `REF_tot` | Total reads assigned to the REF allele |
| `ALT_tot` | Total reads assigned to the ALT allele |
| `MI` | Combined methylation index across both alleles |
| `REF_f` | REF allele methylation fraction (0–1, rounded to 3 dp) |
| `ALT_f` | ALT allele methylation fraction (0–1, rounded to 3 dp) |
| `ref_alt_ratio` | REF/ALT read ratio (balance check; expected ≈ 1.0) |
| `mean_methylation` | Cohort mean methylation for this CpG (from `filter_cpgs`) |
| `sd_methylation` | Cohort SD for this CpG (from `filter_cpgs`) |
| `Category` | CpG variance class: `LOWvar`, `MSDvar`, `SDvar`, or `Mvar` |

---

### 3. `meth_summary_{type}_{sample_id}.xlsx` — Per allele methylation summary

One row per (sample, SNP position, CpG position, allele type) combination.

| Column | Description |
|--------|-------------|
| `sample_id` | Sample identifier |
| `snp_pos` | Genomic position of the phasing SNP |
| `cpg_pos` | Genomic position of the CpG |
| `DMR` | Imprinted DMR name |
| `allele_type` | `"REF"` or `"ALT"` |
| `total_reads` | Total reads for this allele at this CpG |
| `methylated` | Methylated read count |
| `unmethylated` | Unmethylated read count |
| `meth_frac` | Methylation fraction (methylated / total_reads, rounded to 3 dp) |
| `mean_methylation` | Cohort mean methylation for this CpG (from `filter_cpgs`) |
| `sd_methylation` | Cohort SD for this CpG (from `filter_cpgs`) |
| `Category` | CpG variance class: `LOWvar`, `MSDvar`, `SDvar`, or `Mvar` |

---

### 4. `lineplot_{type}_{sample_id}.pdf` — DMR methylation line plots

One page per DMR. Each plot shows `REF_f` and `ALT_f` (REF and ALT allele methylation fractions) across all CpG positions within the DMR, faceted by SNP. Points are shaped by CpG `Category`. Expected pattern for a normally imprinted DMR: one allele near 100% methylation, the other near 0%.

---

### `prepare_cpg_snp_input()` output — `cpg_snps_CG_{type}_{id}.xlsx`

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `pos` | CpG position |
| `context` | Always `"CG"` |
| `total_meth` | Watson + Crick methylated read count |
| `total_cov` | Watson + Crick total coverage |
| `meth_frac` | total_meth / total_cov |
| `DMR` | Imprinted DMR name |
| `snp_pos` | Position of the linked heterozygous SNP |
| `REF` | Reference allele at the SNP |
| `ALT` | Alternative allele at the SNP |
| `GT` | Genotype (always `"0/1"` — heterozygous only) |
| `AD` | Allelic depth string (e.g. `"15,13"`) |
| `DP` | Total SNP read depth |
| `ref_depth` | REF allele read depth |
| `alt_depth` | ALT allele read depth |
| `total_depth` | ref_depth + alt_depth |
| `sample_id` | Sample identifier |

---

## Supported Imprinted DMRs

The package covers the **41 canonical human imprinted DMRs** on GRCh38, including:

| DMR | Chromosome | Associated disorder |
|-----|-----------|---------------------|
| H19/IGF2 | chr11p15.5 | BWS (hypometh) / SRS (hypermeth) |
| KCNQ1OT1 | chr11p15.5 | BWS (hypometh) |
| SNRPN | chr15q11-q13 | PWS / AS |
| MEG3/DLK1 | chr14q32 | Temple syndrome / Kagami-Ogata syndrome |
| PLAGL1 | chr6q24 | Transient neonatal diabetes mellitus |
| GRB10 | chr7p12 | SRS |
| DIRAS3 | chr1p31 | — |
| PPIEL | chr1p36 | — |
| *...and 33 more* | | |

---

## System Requirements

- **R** ≥ 4.1.0
- **Bioconductor** ≥ 3.14
- **samtools** accessible on `PATH` — required at runtime by `extract_bam_regions()` (calls `samtools view`, `samtools sort`, `samtools index`)
- **bgzip / tabix** — **not** required (SNP files are read as plain-text VCF)

---

## Citation

If you use ImprintCapASM in your research, please cite:

> Saha S. *et al.* (2026). ImprintCapASM: SNP-phased allele-specific methylation analysis
> for imprint disorder diagnostics. *R package version 0.1.0*.
> https://CRAN.R-project.org/package=ImprintCapASM

---

## License

MIT © Subham Saha

---

## Contributing

Bug reports and feature requests are welcome via [GitHub Issues](https://github.com/19-saha/ImprintCapASM/issues).
Pull requests should be submitted against the `dev` branch.
