# nf-core/tumourevo: Output

## Introduction

This document describes the output produced by the pipeline. All plots generated in each step are summarised into the final report.
The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Variant Annotation](#variant-annotation)
- [Formatter](#formatter)
- [Lifter](#lifter)
- [Catalogue Driver Annotation](#driver-annotation)
- [QC](#qc)
- [Subclonal Deconvolution](#subclonal-deconvolution)
- [Signature Deconvolution](#signature-deconvolution)

## Directory Structure

The default directory structure is as follows:

```
{outdir}
├── variant_annotation
|   └── vep
│       └── <sample>
├── driver_annotation
|   └── annotate_driver
│       └── <sample>
├── pipeline_info
├── formatter
│   ├── cna2cnaqc
│       └── <sample>
│   ├── cnaqc2tsv
│       └── <patient>
|   └── vcf2cnaqc
│       └── <sample>
├── lifter
│   ├── mpileup
│       └── <patient>
│   └── positions
│       └── <sample>
├── QC
│   ├── tinc
│       └── <sample>
│   ├── CNAqc
│       └── <sample>
│   └── join_CNAqc
│       └── <patient>
├── signature_deconvolution
|   ├── SigProfiler
│       └── <dataset>
|   └── SparseSignatures
│       └── <dataset>
└── subclonal_deconvolution
|   ├── mobster
│       └── <sample>
|   ├── viber
│       └── <patient>
|   ├── ctree
│       └── <patient>,<sample>
|   └── pyclonevi
│       └── <patient>
work/
.nextflow.log
```

<!--Intermediate steps connetting the main subworkflows will output [unpublished results](#unpublished-results) which will be available in the working directory of the pipeline. These steps consist in-->

## Variant Annotation

This directory contains results from the variant annotation subworkflow. At the level of individual samples, genomic variants present in the input VCF files are annotated using [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html). <!--and further converted into MAF format using [vcf2maf](https://github.com/mskcc/vcf2maf). Genomic information of samples from the same cohort are summarized into a unique MAF object with [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html).-->

### VEP

[VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html) is a `Ensembl` tool that determines the effect of all variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence.
This step starts from VCF files.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/variant_annotation/vep/<dataset>/<patient>/<sample>/`**

- `<dataset>_<patient>_<sample>.vcf.gz` and `<dataset>_<patient>_<sample>.vcf.gz.tbi`
  - VCF file and tabix index with called mutations

</details>

<!-- ### vcf2maf

[vcf2maf](https://github.com/mskcc/vcf2maf) convert a VCF file into a MAF (Mutation Annotation Format) file, where each variant is mapped to only one of all possible gene transcripts/isoforms that it might affect. vcf2maf is designed to work with VEP annotation output.

> **NB:** While VEP is tolerant of chromosome format mismatches (when the input .VCF file uses the UCSC format chrN and the reference fasta uses Ensembl/NCBI format N), vcf2maf is not.
> Make sure the reference fasta chromosome format matches that of your input

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/VariantAnnotation/VCF2MAF/<dataset>/<patient>/<sample>/`**
* `data_vep.maf`
    * annotated MAF file

</details>

### maftools

[maftools]( https://bioconductor.org/packages/release/bioc/html/maftools.html) is an R-based tool that provides a comprehensive set of functions for processing MAF files and to perform most commonly used analyses in cancer genomics. In particular, maftools summarize, analyze and visualize MAF files from the same cohort providing summary statistics of the top mutated genes and aggregating metadata into a single oncoplot visualization.
This step started from annotated MAF files.

MAF fields requirements:

* Mandatory fields: `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`, `Variant_Type` and `Tumor_Sample_Barcode`.

* Optional fields: `VAF` (Variant Allele Frequency), `amino acid change` information.

<details markdown="1">
<summary>Output files for the dataset</summary>

**Output directory: `{outdir}/VariantAnnotation/MAFTOOLS/<dataset>/`**
* `maf_merged.rds`
    * summarized MAF object
* `maf_summary.pdf`
    * summary statistics plot
* `oncoplot.pdf`
    * oncoprint plot

</details> -->

## Formatter

The Formatter subworkflow is used to convert file to other formats and to standardize the output files resulting from different mutation (Mutect2, Strelka) and cna callers (ASCAT,Sequenza).

### vcf2cnaqc

This parser is designed to process VCF files generated by various variant calling tools and to convert them into a unified RDS file format.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/formatter/vcf2cnaqc/<dataset>/<patient>/<sample>/`**

- `<dataset>_<patient>_<sample>_snv.rds`
  - RDS file containing parsed VCF in table format

</details>

### cna2cnaqc

This parser is designed to standardize copy number calls and purity estimates from various callers into a unified format.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/formatter/cna2cnaqc/<dataset>/<patient>/<sample>/`**

- `<dataset>_<patient>_<sample>_cna.rds`
  - RDS file containing parsed segments and purity estimate output in table format

</details>

### cnaqc2tsv

This parser is designed to convert mutations data of joint CNAqc analysis from CNAqc format (RDS file) into a tabular format (TSV file). This step is mandatory for running python-based tools (e.g. PyClone-VI, SigProfiler) and it is mandatory if `--tools` contains either `pyclone-vi` or `sigprofiler`.

<details markdown="1">
<summary>Output files for all patients</summary>

**Output directory: `{outdir}/formatter/cnaqc2tsv/<dataset>/<patient>/`**

- `<dataset>_<patient>_joint_table.tsv`
  - TSV file containing mutations mapped to corrsponding copy number segments.

</details>

## Lifter

The Lifter subworkflow is an optional step and it is run when single sample VCF file are provided. When multiple samples from the same patient are provided, the user can specify either a single joint VCF file, containing variant calls from all tumor samples of the patient (see [joint variant calling](<(https://nf-co.re/sarek/3.4.2/parameters/#joint_mutect2)>)), or individual sample specific VCF files. In the latter case, path to tumor BAM files must be provided in order to collect all mutations from the samples and perform pile-up of sample's private mutations in all the other samples. Two intermediate steps, [get_positions](#get_positions) and [mpileup](#mpileup), are performed to identify private mutations in all the samples and retrieve their variant allele frequency. Once private mutations are properly defined, they are merged back into the original VCF file during the [join_positions](#join_positions) step. The updated VCF file is then converted into a RDS object.

### mpileup

At this stage, [bcftools](https://samtools.github.io/bcftools/bcftools.html) mpileup is run to retrieve frequency information of private mutations across all samples.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/lifter/mpileup/<dataset>/<patient>/<sample>/`**

- `<dataset>_<patient>_<sample>.bcftools_stats.txt`
  - TXT file with statistics on the called mutations
- `<dataset>_<patient>_<sample>.vcf.gz` and `<dataset>_<patient>_<sample>.vcf.gz.tbi`:
  - VCF file and tabix index with called mutations

</details>

### positions

This step allows to retrieve private and shared mutations across samples originated from the same patient. Previously retrieved mutations are joined with original mutations present in input VCF, which is in turn converted into an RDS object using [vcfR](https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html).

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/lifter/positions/<dataset>/<patient>/<sample>/`**

- `<dataset>_<patient>_<sample>.pileup_VCF.rds`
  - RDS containing shared and private mutations
- `<dataset>_<patient>_<sample>.positions_missing`
  - TXT file containing mutations to be retrieved for a given sample

</details>

<details markdown="1">
<summary>Output files for all patients</summary>

**Output directory: `{outdir}/lifter/positions/<dataset>/<patient>/`**

- `<dataset>_<patient>__all_positions.rds`
  - RDS containing shared and private mutations

</details>

## Catalogue Driver Annotation

This directory contains results from the driver annotation subworkflow.

### Tumour-type driver annotation

According to the specified tumour type, potential driver mutations are identified and annotated using [IntOGen database](<(https://www.intogen.org/search)>). The user can also provide a custom table of driver genes.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/driver_annotation/annotate_driver/<dataset>/<patient>/<sample>/`**

- `<dataset>_<patient>_<sample>_driver.rds`
  - RDS with annotated mutations

</details>

<!-- ### Build reference

Add description

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/driver_annotation/BuildReference/<dataset>/<patient>/<sample>/`**
* `fit.rds`
    * add description
* `plot.rds`
    * add description

</details>

### dndsCV

Add description

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/driver_annotation/DNDSCV/<dataset>/<patient>/<sample>/`**
* `dnds.rds`
    * DNDSCV object in RDS format

</details> -->

## QC

The QC subworkflows requires in input a segmentation file from allele-specific copy number callers (either [Sequenza](https://sequenzatools.bitbucket.io/#/home), [ASCAT](https://github.com/VanLoo-lab/ascat)) and the joint VCF file. As a first step, the QC subworkflow provides an estimate of normal and tumour samples contamination in [TINC](#tinc) step, in order to have a measure of experimental quality. Then,it first conducts a quality control on copy number and somatic mutation data for individual samples in [CNAqc](#cnaqc) step, and subsequently summarize validated information at patient level in [join_CNAqc](#join_cnaqc) step.
The QC subworkflow is a crucial step of the pipeline as it ensures high confidence in identifying clonal and subclonal events while accounting for variations in tumor purity.

### TINC

[TINC](https://caravagnalab.github.io/TINC/index.html) is a package to calculate the contamination of tumor DNA in a matched normal sample. TINC provides estimates of the proportion of cancer cells containing the normal sample (TIN), and the proportion of cancer cells in the tumor sample (TIT).

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/QC/tinc/<dataset>/<patient>/<sample>/`**

- `<dataset>_<patient>_<sample>_fit.rds`
  - TINC fit containing TIN and TIT estimates in RDS;
- `<dataset>_<patient>_<sample>_plot.rds` and `<dataset>_<patient>_<sample>_plot.pdf`:
  - TINC report contianign TIN and TIT plots in PDF and RDS;
- `<dataset>_<patient>_<sample>_qc.csv`:
  - TINC summary report on normal contamination.

</details>

### CNAqc

[CNAqc](https://caravagnalab.github.io/CNAqc/) is a package that performs quality control of bulk cancer sequencing data for validating copy number segmentations against variant allele frequencies of somatic mutations.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/QC/CNAqc/<dataset>/<patient>/<sample>/`**

- `<dataset>_<patient>_<sample>_data_plot.rds` and `<dataset>_<patient>_<sample>_data.pdf`
  - CNAqc report with genome wide mutation and allele specific copy number plots in RDS and PDF
- `<dataset>_<patient>_<sample>_qc_plot.rds` and `<dataset>_<patient>_<sample>_qc.pdf`
  - QC step report resulting from peak analysis in RDS and PDF
- `<dataset>_<patient>_<sample>_qc.rds`
  - CNAqc RDS object

</details>

### join_CNAqc

This module creates a multi-CNAqc object for patient by summarizing the quality check performed at the single sample level. For more information about the structure of multi-CNAqc object see [CNAqc documentation](<(https://caravagnalab.github.io/CNAqc/)>).

<details markdown="1">
<summary>Output files for all patients</summary>

**Output directory: `{outdir}/QC/join_CNAqc/<dataset>/<patient>/`**

- `<dataset>_<patient>_multi_cnaqc_ALL.rds`
  - unfiltered mCNAqc RDS object
- `<dataset>_<patient>_multi_cnaqc_PASS.rds`
  - filtered mCNAqc RDS object

</details>

## Signature Deconvolution

<!-- Mutational signatures represent characteristic patterns of somatic mutations in cancer genomes, reflecting the underlying mutational processes at the basis of tumor evolution and progression. Mutational signatures are discovered by analyzing ensemble point-mutation counts from a set of individual samples. Validated mutations from [join_CNAqc](#join_cnaqc) step are converted into a TSV joint table in (see [tsvparse](#tsvparse) module), subsequently given as input to signature deconvolution subworkflow, which performs de novo extraction, inference, deciphering or deconvolution of mutational counts.  -->

Mutational signatures are distinctive patterns of somatic mutations in cancer genomes that reveal the underlying mutational processes driving tumor evolution and progression. These signatures are identified by analyzing aggregated point-mutation counts from multiple samples. Validated mutations from the [join_CNAqc](#join_cnaqc) step are converted into a joint TSV table (see [cnaqc2tsv](#cnaqc2tsv)) and then input into the signature deconvolution subworkflow, which performs de novo extraction, inference, interpretation, or deconvolution of mutational counts.

The results of this step are collected in `{pubslish_dir}/signature_deconvolution/`. Two tools can be specified by using `--tools` parameter: [SparseSignatures](#sparsesignatures) and [SigProfiler](#sigprofiler).

### SparseSignatures

[SparseSignatures](https://www.bioconductor.org/packages/release/bioc/html/SparseSignatures.html) is an R-based computational framework that provides a set of functions to extract and visualize the mutational signatures that best explain the mutation counts of a large number of patients.

<details markdown="1">
<summary>Output files for dataset</summary>

**Output directory: `{outdir}/signatures_deconvolution/SparseSig/<dataset>/`**

- `<dataset>_best_params_config.rds`
  - signatures best configiration object
- `<dataset>_cv_means_mse.rds`
  - cross validation output RDS
- `<dataset>_nmf_Lasso_out.rds`
  - NMF Lasso output RDS
- `<dataset>_plot_signatures.pdf` and `<dataset>_plot_signatures.rds`
  - exposure plot in PDF and RDS

</details>

### SigProfiler

[SigProfiler](https://osf.io/t6j7u/wiki/home/) is a python framework that allows _de novo_ extraction of mutational signatures from data generated in a matrix format. The tool identifies the number of operative mutational signatures, their activities in each sample, and the probability for each signature to cause a specific mutation type in a cancer sample. The tool makes use of `SigProfilerMatrixGenerator` and `SigProfilerExtractor`, seamlessly integrating with other `SigProfiler` tools.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/signatures_deconvolution/SigProfiler/<dataset>/results`**

- `input/`
  - folder containing a copy of the user-provided input files for SigProfilerMatrixGenerator step
- `input_data.txt`
  - join table of all mutations in the dataset in TXT
- `logs/`
  - folder containing the error and log files for SigProfilerMatrixGenerator step
- `output/`
  - folder containing the DBS, SBS, INDEL nucleotide matrices resulting from SigProfilerMatrixGenerator step
- `SBS96/`
  - folder containing the results of SigProfilerExtractor step in the SBS96 mutational context. This directory will contain:
    - `All_Solutions/`
      - subdirectory containing the results from running extractions at each rank within the range of the input. For more details visit the [official website](https://osf.io/t6j7u/wiki/5.%20Output%20-%20All%20Solutions/)
    - `Suggested_Solution/`
      - subdirectory containing the optimal solution. For more details visit the [official website](https://osf.io/t6j7u/wiki/6.%20Output%20-%20Suggested%20Solution/)
- `JOB_METADATA.txt`
  - TXT file containing all the metadata about the system and runtime of the job
- `Seeds.txt`
  - TXT file containing the replicate IDs and preset seeds

</details>

## Subclonal Deconvolution

The subclonal deconvolution subworkflow requires in input a joint `mCNAqc` object resulting from the [join_CNAqc](#join_cnaqc) step. The subworkflow will perform multi-sample deconvolution if more than one sample for each patient is present.

<!-- The user can run MOBSTER before performing subclonal deconvolution in three different ways by specifing the parameter `--remove_tail`. -->

<!-- > **NB:** If `--remove_tail all,once` and MOBSTER is not specified in the tools for subclonal deconvolution, the pipeline will throw an error. -->

<!-- Two different modalities can be specified by the user using `--mode` parameter: single sample and multi sample mode. -->
<!-- If `--mode singlesample` is provided, each sample is analysed individually providing a snapshot of clonal and subclonal diversity starting from allele frequency of detected somatic variants. When multiple samples from the same patient are provided, the user may take advantage of the multisample modality, by setting `--mode multisample`. -->
<!-- This approach allows for a more detailed and accurate identification of subclonal populations, as it can capture spatial and temporal heterogeneity within the tumor. By integrating data across multiple samples, it improves the resolution of subclonal structures and provides insights into the evolutionary dynamics and progression of the tumor. -->

<!-- Various tools can be specified using the `--tools` parameter, leading to different methods for performing subclonal deconvolution analysis. Among the available tools, [MOBSTER](https://caravagnalab.github.io/mobster/) operates only on single samples but can still be specified for use in multi-sample mode, while [VIBER](https://caravagnalab.github.io/VIBER/index.html) and [PyClone-VI](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03919-2) can operates in both modalities. In case `--mode multisample` and `--tools mobster,pyclone-vi` are specified, first MOBSTER is run on individual samples to remove tail mutations and then PyClone-VI operates a multivariate subclonal deconvolution on the preprocessed MOBSTER mutations. A similar procedure is perfomed when  `--mode multisample` and `--tools mobster,viber` are defined. More detailed explanation is provided in the following sections.  -->

The results of subclonal decovnultion step are collected in `{outdir}/subclonal_deconvolution/` directory.

<!-- ### Single sample

If `--mode singlesample` is provided, each sample is analysed individually providing a snapshot of clonal and subclonal diversity starting from allele frequency of detected somatic variants. Available tools provides different ways of performing subclonal deconvolution. Both MOBSTER and VIBER model mutation counts as mixture of binomial distribution; however, MOBSTER includes a pareto Type-I power law distribution to model within-clone neutral dynamics. Instead, PyClone-VI models clonal architecture taking into account variant allele frequencies corrected for coincident copy number variation. -->

### MOBSTER

[MOBSTER](https://caravagnalab.github.io/mobster/) is a package that models mutant allelic frequencies and copy-number status by integrating evolutionary theory and Bayesian proabilistic modelling to identify clusters of variants with similar cellular proportions. Futhermore, MOBSTER models the dynamics of passenger mutations via a Pareto distribution giving rise to the so called neutral tail.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/mobster/<dataset>/<patient>/<sample>/`**

- `<dataset>_<patient>_<sample>_mobsterh_st_fit.rds`
  - RDS object contains all fits of subclonal deconvolution
- `<dataset>_<patient>_<sample>_mobsterh_st_best_fit.rds`
  - RDS object contains best fit of subclonal deconvolution
- `<dataset>_<patient>_<sample>_mobsterh_st_best_fit_plots.rds`
  - summary plots of best fits in RDS
- `<dataset>_<patient>_<sample>_REPORT_plots_mobster.{rds,png,pdf}`
  - report of mobster deconvolution in RDS, PDF and PNG format

</details>

### PyClone-VI

[PyClone-VI](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03919-2) is a computationally efficient Bayesian statistical method for inferring the clonal population structure of cancers, by considering allele fractions and coincident copy number variation using a variational inference approach. It works for patients with both single and multiple samples.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/pyclonevi/<dataset>/<patient>`**

- `<dataset>_<patient>_pyclone_input.tsv`
  - TSV file with Pyclone-VI input table
- `<dataset>_<patient>_all_fits.h5`
  - HDF5 file for all possible fit and summary stats
- `<dataset>_<patient>_best_fit.txt`
  - TSV file for the best fit
- `<dataset>_<patient>_cluster_table.csv`
  - CSV file with clone assignment

</details>

### VIBER

[VIBER](https://caravagnalab.github.io/VIBER/index.html) is an R package that implements a variational Bayesian model to fit multi-variate Binomial mixtures. In the context of subclonal deconvolution VIBER models read counts that are associated with the most represented karyotype. It works for patients with both single and multiple samples.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/viber/<dataset>/<patient>`**

- `<dataset>_<patient>_viber_best_st_fit.rds`
  - RDS file for best standard fit
- `<dataset>_<patient>_viber_best_st_fit_plots.rds`
  - RDS file containing summary plots for best standard fit
- `<dataset>_<patient>_viber_best_st_heuristic_fit.rds`
  - RDS file for best standard fit with applied heuristic
- `<dataset>_<patient>_viber_best_st_heuristic_fit_plots.rds`
  - RDS file containing summary plots for best standard fit with applied heuristic
- `<dataset>_<patient>_viber_best_st_mixing_plots.rds`
  - RDS file containing mixing proportion plot for best standard fit
- `<dataset>_<patient>_viber_best_st_heuristic_mixing_plots.rds`
  - RDS file containing mixing proportion plot for best standard fit with applied heuristic
- `<dataset>_<patient>_REPORT_plots_viber.{rds,pdf,png}`
  - report of VIBER deconvolution in RDS, PDF and PNG

</details>

<!-- ### Multi samples

When multiple samples for the same patient are available (e.g. multi-regional or longitudinal experiment), we can be intersted in visualizing the tumor evolution at different level. In multi sample mode, a patient-specific subclonal deconvolution is perfomed.

> **NB:** At this stage, it is strongly recommended that the user conducts single-sample subclonal deconvolution using `mosbter` to eliminate mutations assigned to the neutral tail before proceeding with multi-variate clone identification.
> This process ensures a cleaned signal for downstream analyses aimed at focusing on functional intratumor heterogeneity.  -->

<!-- #### Pyclone-VI

This folder contains the results of multivariate analysis using Pyclone-VI, which can be run prior to removal of mutations assigned to tail in all the samples.

<details markdown="1">
<summary>Output files for all patients with MOBSTER</summary>

**Output directory: `{outdir}/subclonal_deconvolution/pyclonevi/<dataset>/<patient>/`**

- `<patient>_with_mobster_all_fits.h5`
  - HDF5 file for all possible fit
- `<patient>_with_mobster_best_fit.txt`
  - TSV file for the best fit

- `all_fits.h5`
  - HDF5 file for all possible fit and summary stats
- `best_fit.txt`
  - TSV file for the best fit

</details>

<details markdown="1">
<summary>Output files for all patients without MOBSTER</summary>

**Output directory: `{outdir}/subclonal_deconvolution/pyclonevi/<dataset>/<patient>/`**

- `without_mobster_all_fits.h5`
  - HDF5 file for all possible fit and summary stats
- `without_mobster_best_fit.txt`
  - TSV file for the best fit

</details> -->

<!-- #### VIBER

This folder contains the results of multivariate analysis using VIBER, which can be run prior to removal of mutations assigned to tail in all the samples.

<details markdown="1">
<summary>Output files for all patients with MOBSTER</summary>

**Output directory: `{outdir}/subclonal_deconvolution/viber/<dataset>/<patient>/`**

- `viber_best_st_fit.rds`
  - RDS file for best standard fit
- `viber_best_st_mixing_plots.rds`
  - RDS file containing mixing proportion for best standard fit
- `viber_best_st_heuristic_fit.rds`
  - RDS file for best standard fit with applied heuristic
- `viber_best_st_heuristic_mixing_plots.rds`
  - RDS file containing mixing proportion for best standard fit  with applied heuristic

</details>

<details markdown="1">
<summary>Output files for all patients without MOBSTER</summary>

**Output directory: `{outdir}/subclonal_deconvolution/viber/<dataset>/<patient>/`**

- `viber_without_mobster_best_st_fit.rds`
  - RDS file for best standard fit
- `viber_without_mobster_best_st_mixing_plots.rds`
  - RDS file containing mixing proportion for best standard fit
- `viber_without_mobster_best_st_heuristic_fit.rds`
  - RDS file for best standard fit with applied heuristic
- `viber_without_mobster_best_st_heuristic_mixing_plots.rds`
  - RDS file containing mixing proportion for best standard fit  with applied heuristic

</details> -->

### ctree

Subclonal deconvolution results are used to build clone tree from both single samples and multple samples using [ctree](https://caravagnalab.github.io/ctree/index.html). ctree is a R-based package which implements basic functions to create, manipulate and visualize clone trees by modelling Cancer Cell Fractions (CCF) clusters. Annotated driver genes must be provided in the input data.

> **NB:** When `--tools pyclone-vi` is used, the output of PyClone-VI subclonal deconvolution is preprocessed prior to clone tree inference. Since ctree requires labeling one of the clusters as "clonal," the one with the highest CCF across all samples is choosen.

VIBER and MOBSTER fits are already compatible for ctree analysis.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/subclonal_deconvolution/ctree/<dataset>/{<patient>,<patient>/<sample>/}`**

- `{<dataset>_<patient>,<dataset>_<patient>_<sample>}_ctree_<tool>.rds`
  - RDS file containing inferred clone tree
- `{<dataset>_<patient>,<dataset>_<patient>_<sample>}_ctree_<tool>_plots.rds`
  - RDS file for clone tree plot
- `{<dataset>_<patient>,<dataset>_<patient>_<sample>}_REPORT_plots_ctree_<tool>.{rds,png,pdf}`
  - ctree report in RDS,PNG and PDF

</details>
<!--
### Multi sample
<details markdown="1">
<summary>Output files for all patients</summary>

**Output directory: `{outdir}/subclonal_deconvolution/ctree/<patient>/`**

- `ctree_<tool>.rds`
  - RDS file containing inferred clone tree
- `ctree_<tool>_plots.rds`
  - RDS file for single sample clone tree plot
- `ctree_input_pyclonevi.csv`
  - CSV file required for clone tree inference from pyclone

</details> -->

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

## Reference files

Different tools of the pipeline generate references files. Once reference file for VEP and SigProfiler are not provided they are stored in the tool-specific folder.

### VEP

When VEP cache is not specified, the desired VEP cache is downladed in `{outdir}/references/VEP/vep_cache/homo_sapiens/{VEP_version}_{ref_genome}`.

### SigProfiler

Reference genome for SigProfiler is store in the following folder:

<details markdown="1">
<summary>Output files</summary>

**Output directory: `{outdir}/subclonal_deconvolution/signature_deconvolution/SigProfiler/genome/tsb/{ref_genome}`**

- `{chromosome}.txt`
  - genome assemby chromosme level
- `{ref_genome}_proportions.txt`
  - genome assemby proportions
  </details>
