# nf-core/tumourevo: Usage

## Table of contents

- [nf-core/tumourevo: Usage](#nf-coretumourevo-usage)
    - [Table of contents](#table-of-contents)
- [Introduction](#introduction)
- [Running the pipeline](#running-the-pipeline)
    - [Quickstart](#quickstart)
    - [Input: Sample sheet configurations](#input-sample-sheet-configurations)
        - [Overview: Samplesheet Columns](#overview-samplesheet-columns)
    - [Pipeline modalities](#pipeline-modalities)
        - [Variant calling](#variant-calling)
            - [Multi-sample variant calling](#1-multi-sample-variant-calling)
            - [Single sample variant calling](#2-single-sample-variant-calling)
        - [Subclonal deconvolution](#subclonal-deconvolution)
    - [Driver annotation](#driver-annotation)
    - [Avaiable tools](#avaiable-tools)
        - [Updating the pipeline](#updating-the-pipeline)
        - [Reproducibility](#reproducibility)
    - [Main arguments](#main-arguments)
        - [`-profile`](#-profile)
    - [Job resources](#job-resources)
        - [Automatic resubmission](#automatic-resubmission)
        - [Custom resource requests](#custom-resource-requests)
    - [Other command line parameters](#other-command-line-parameters)
        - [`--outdir`](#--outdir)
        - [`--email`](#--email)
        - [`--email_on_fail`](#--email_on_fail)
        - [`--max_multiqc_email_size`](#--max_multiqc_email_size)
        - [`-name`](#-name)
        - [`-resume`](#-resume)
        - [`-c`](#-c)
        - [`--custom_config_version`](#--custom_config_version)
        - [`--custom_config_base`](#--custom_config_base)
        - [`--max_memory`](#--max_memory)
        - [`--max_time`](#--max_time)
        - [`--max_cpus`](#--max_cpus)
        - [`--plaintext_email`](#--plaintext_email)
        - [`--monochrome_logs`](#--monochrome_logs)
        - [`--multiqc_config`](#--multiqc_config)

# Introduction

**tumourevo** is a workflow to infer a tumour evolution model from whole-genome sequencing (WGS) data.

Through the analysis of variant and copy-number calls, it reconstructs the evolutionary process leading to the observed tumour genome. Most of the analyses can be done at mutliple levels: single sample, multiple samples from the same patient (multi-region/longitudinal assays), and multiple patients from distinct cohorts.

# Running the pipeline

## Quickstart

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/tumourevo \
 -r <VERSION> \
 -profile <PROFILE> \
 --input <INPUT CSV> \
 --outdir ./results
 --tools <TOOLS>
```

`-r <VERSION>` is optional but strongly recommended for reproducibility and should match the latest version.

`-profile <PROFILE>` is mandatory and should reflect either your own institutional profile or any pipeline profile specified in the [profile section](##-profile).

This documentation imply that any `nextflow run nf-core/tumourevo` command is run with the appropriate `-r` and `-profile` commands.

This will launch the pipeline and perform variant calling with the tools specified in `--tools`, see the [parameter section]([https://github.com/caravagnalab/nf-core-tumourevo/tree/dev]) for details on the available tools.

Unless running with the `test` profile, the paths of input files must be provided within the `<INPUT CSV>` file specified in `--input`, see the [input section]([https://github.com/caravagnalab/nf-core-tumourevo/tree/dev]) for input requirements.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Input: Sample sheet configurations

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the parameter `--input` to specify its location. It has to be a comma-separated file with at least 5 columns, and a header row as shown in the examples below.

It is recommended to use the absolute path of the files, but a relative path should also work.

For the joint analysis of multiple samples, a tumor BAM file is required for each sample, such that the number of reads of a private mutation can be retrieved for all the samples thorugh `mpileup`.

Multiple samples from the same patient must be specified with the same `dataset` ID, `patient` ID, and a different `tumour_sample` ID. `normal_sample` columns is required.

Multiple patients from the same dataset must be specified with the same `dataset` ID, and a different `patient` ID.

**tumourevo** will output sample-specific results in a different directory for _each sample_, patient-specific results in a common directory for _each patient_, and dataset-specific results in a common directory for _each dataset_.

Output from different workflows, subworkflows and modules will be in a specific directory for each dataset, patient, sample and tool configuration.

### Overview: Samplesheet Columns

| Column    | Description                                                                                                                                                                                                                                                                                                                       |
| --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dataset` | **Dataset ID**; when sequencing data from multiple datasets is analysed, it designates the source dataset of each patient; must be unique for each dataset, but one dataset can contain samples from multiple patients. <br /> _Required_                                                                                      |
| `patient` | **Patient ID**; designates the patient/subject; must be unique for each patient, but one patient can have multiple samples (e.g. from multiple regions or multiple time points). <br /> _Required_                                                                                                                                |
| `tumour_sample`  | **Sample ID** for each sample; more than one sample for each subject is possible. Must match the sample ID present in the VCF.  <br /> _Required_
| `normal_sample` | **Normal sample ID** of each sample. Must match the normal sample ID present in the VCF. <br /> _Required_      |                                 |
| `vcf`  | Full path to the vcf file. <br /> _Required_                                                                                                        |
| `tbi`  | Full path to the vcf `tabix` index file. <br /> _Required_                                                                                      |
|`cna_caller`| Name of the copy number caller used to generate your data. <br /> _Required_ |
| `cna_segments` | Full path to the segmentation files and copy number state from copy-number calling. <br /> _Required_ |
| `cna_extra` | Full path to files including the ploidy and purity estimate from the copy-number caller. <br /> _Required_ |
| `cancer_type` | Tumour type (either `PANCANCER` or one of the tumor type present in the driver table) <br /> _Required_ |
| `tumour_alignment`  | Full path to the tumour bam file. <br /> _Optional_                                                       |
| `tumour_alignment_index`  | Full path to the tumour bam index file. <br /> _Optional_                                              |

An [example samplesheet](https://github.com/caravagnalab/tumourevo/blob/dev/test_input.csv) has been provided with the pipeline.

## Pipeline modalities

The tumourevo pipeline supports variant annotation, driver annotation, quality control processes, subclonal deconvolution and signature deconvolution analysis through various tools. It can be used to analyse both single sample experiments and longitudinal/multi-region assays, in which multiple samples of the same patient are avaiable.
As input, you must provide at least information on the samples, the VCF file from one of the supported callers and the output of one of the supported copy number caller. By default, if multiple samples from the same patient are provided, they will be analysed in a multivariate framework (which affects in particular the subclonal deconvolution deconvolution steps) to retrieve information useful in the reconstruction of the evolutionary process. Depending on the variant calling strategy (single sample or multi sample) and the provided input files, different strategies will be applied.

<!-- aggiungi un riassunto di cosa voglia dire single e multi sample (analisi multivariata, soprattutto per subclonal deconv)
E' possibile usarla sia nel caso di vc multi sample che indipendente -->

## Variant calling

### 1. Multi-sample variant calling

Modern tools (ie: Platypus and Mutect2) allow to perform variant calling directly in multisample mode. If the VCFs provided as input are already multisample, no additional step is required.

<!-- If the variant calling had been performed indepentently on each sample from the same patient, 

If you run the pipeline in `singlesample` mode, all the samples, even if belonging to the same patient, are assumed to be independent. In this framework, the subclonal deconvolution is affected, identifying clonal and subclonal composition of the sample starting from allele frequency of detected somatic variants and identifying mutagenic processes for each independent element.  -->

<!-- rephrase better -->
<!-- Si assume che i campioni siano indipendenti e che quindi la subclonal deconv viene svolta cercando popolazioni sottoclonali a lv di singolo campione. Sign deconv si cercano i processi mutagenici comuni in un dataset fatto di elementi indipenti -->

#### Examples

Running the pipeline

```bash
nextflow run nf-core/tumourevo \
 -r <VERSION> \
 -profile <PROFILE> \
 --input <INPUT CSV> \
 --outdir results \
 --tools pyclonevi,mobster,viber,sparsesignature,sigprofiler
```

Minimal input file, two samples from the same patient:

```bash
dataset,patient,sample,normal_sample,vcf,vcf_tbi,cna_segments,cna_extra,cna_caller,cancer_type
dataset1,patient1,S1,N1,patient1_S1.vcf.gz,patient1_S1.vcf.gz.tbi,/CNA/patient1/S1/segments.txt,/CNA/patient1/S1/purity_ploidy.txt,caller,PANCANCER
dataset1,patient1,S2,N1,patient1_S1.vcf.gz,patient1_S1.vcf.gz.tbi,/CNA/patient1/S2/segments.txt,/CNA/patient1/S2/purity_ploidy.txt,caller,PANCANCER
```

### 2. Single sample variant calling

If the variant calling performed independently on each sample, even if coming from the same patient, you can provide the BAM and BAI files from each tumor sample. In this way, a classical pileup strategy will be used in order to retrieve the depth for all samples of private mutations, in order to correctly perform the subclonal deconvolution analysis.

Input file for two patients without joint variant calling, bam files available:

```bash
dataset,patient,sample,normal_sample,vcf,vcf_tbi,cna_segments,cna_extra,cna_caller,cancer_type,tumour_bam,tumour_bai
dataset1,patient1,S1,N1,patient1_S1.vcf.gz,patient1_S1.vcf.gz.tbi,/CNA/patient1/S1/segments.txt,/CNA/patient1/S1/purity_ploidy.txt,caller,PANCANCER,patient1/BAM/S1.bam,patient1/BAM/S1.bam.bai 
dataset1,patient1,S2,N1,patient1_S2.vcf.gz,patient1_S2.vcf.gz.tbi,/CNA/patient1/S2/segments.txt,/CNA/patient1/S2/purity_ploidy.txt,caller,PANCANCER,,patient1/BAM/S2.bam,patient1/BAM/S2.bam.bai  
```

If you can not include the bam files in the input csv, the pipeline will run anyway, treating each sample as independent.

<!-- You can use the `multisample` mode of tumourevo to analyse samples from multi-region and longitudinal assays. This allows you to track in space and time the existing tumor populations, and better understand its heterogeneity. This modality integrates data across multiple samples, thus improving the resolution of subclonal structures and providing insights into the evolutionary dynamics and progression of the tumor.
Two of the avaiable tools for subclonal deconvolution, `pyclonevi` and `viber` can by-design be run in multi-sample mode, inferring the subclonal structure of samples. If you add `mobster` to the `--tool` parameter when running the pipeline in this modality, it will be run at first on each individual sample (since the tool does not support at the moment multi-sample analysis) in order to recognize neutral tail mutations and remove them. The mutations data manipulated in this way will then be processed by either `pyclone`, `viber` or both using the multivariate subclonal deconvolution as described before.  -->

# Subclonal deconvolution

### The `remove_tail` parameter

#### Example

Running the pipeline without mobster

```bash
nextflow run nf-core/tumourevo \
 -r <VERSION> \
 -profile <PROFILE> \
 --samples <INPUT CSV> \
 --publish_dir ./results \
 --mode multisample \
 --tools pyclone,viber,ctree,sparsesignature,sigprofiler
```

Running the pipeline with mobster

<!-- mobster lavora single sample, mentre pyclone e viber anche multisample -->

```bash
nextflow run nf-core/tumourevo \
 -r <VERSION> \
 -profile <PROFILE> \
 --samples <INPUT CSV> \
 --publish_dir ./results \
 --mode multisample \
 --tools mobster,pyclone,viber,ctree,sparsesignature,sigprofiler
```

Minimal input file (two samples from the same patient):

```bash
dataset,patient,sample,vcf,vcf_tbi,cna_dir,cna_caller
dataset1,patient1,S1,patient1_S1.vcf.gz,patient1_S1.vcf.gz.tbi,/CNA/patient1/S1,caller
dataset1,patient1,S2,patient1_S2.vcf.gz,patient1_S2.vcf.gz.tbi,/CNA/patient1/S2,caller
```

## Driver annotation

You can retrieve tumor-specific drivers in the driver annotation step by specifying the tumor type in the input csv. Pan-cancer drivers will be retrieved by specifying `PANCANCER` as tumour type in the input csv file.
For this step, we currently refer to [IntOGen latest release](https://www.nature.com/articles/s41568-020-0290-x), but you can also provide a custom driver table that will be used in the analysis.
Please note that the tumor types reported in the input file must correspond to those present in the table used for the annotation (default driver table used can be found [here](https://github.com/caravagnalab/nextflow_modules/blob/main/2023-05-31_IntOGen-Drivers/Unfiltered_drivers.tsv))

<!-- add a snapshot of the driver table -->

Input file:

```bash
dataset,patient,sample,vcf,vcf_tbi,cna_dir,cna_caller,tumour_bam,tumour_bai,tumour_type
dataset1,patient1,S1,patient1_S1.vcf.gz,patient1_S1.vcf.gz.tbi,/CNA/patient1/S1,caller,patient1/BAM/S1,patient1/BAM/S1.bam.bai,BRCA
dataset1,patient1,S2,patient1_S2.vcf.gz,patient1_S2.vcf.gz.tbi,/CNA/patient1/S2,caller,patient1/BAM/S2,patient1/BAM/S2.bam.bai,BRCA 
dataset1,patient1,S1,patient1_S1.vcf.gz,patient1_S1.vcf.gz.tbi,/CNA/patient2/S1,caller,patient1/BAM/S1,patient1/BAM/S1.bam.bai,BRCA 
dataset1,patient2,S2,patient2_S2.vcf.gz,patient2_S2.vcf.gz.tbi,/CNA/patient2/S2,caller,patient2/BAM/S2,patient2/BAM/S2.bam.bai,BRCA 
```

## Avaiable tools

We report the different tools included in the pipeline.

1. **Gene annotation**

- [EnsemblVEP]()

2. **Driver annotation**
3. **Quality control**

- [CNAqc](https://caravagnalab.github.io/CNAqc/)

<!-- #### TINC

#### INCOMMON -->

4. **Subclonal deconvolution**

- [MOBSTER](https://caravagnalab.github.io/mobster/)
- [PyClone-VI](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03919-2)
- [VIBER](https://caravagnalab.github.io/VIBER/index.html)
- [Ctree](https://caravagnalab.github.io/ctree/)

5. **Signature deconvolution**

- [SparseSignature]()
- [SigProfiler]()
<!-- ### Start with annotation

Example of the values used for these parameters:

```bash
ref_genome_vep =  "$HOME/ref_genomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
vep_dir_cache  =  "$HOME/ref_genomes/VEP/"
vep_cache_version = "110"
vep_species = "homo_sapience" 
assembly = "GRCh38"
```

Using VEP plugins:
`--plugin` - plugin modules should be installed in the Plugin subdirectory of the VEP cache directory (defaults to "$HOME/.vep/").
To enable specify `--plugin SingleLetterAA`.

### Start with Maftools

`read.maf` function reads multiple MAF files (e.g. multisample/multipatient cohort), summarizes it in various ways and stores it as an MAF object.
MAF object contains main maf file, summarized data and any associated sample annotations.
For visualization it uses `plotmafSummary` function to plot the summary of the maf file, which displays number of variants in each sample and variant types summarized by Variant_Classification. Better representation of maf file can be shown as oncoplots using `oncoplot` function. -->

<!-- ### Subclonal Deconvolution

This step can be started either from XXX files or XXXX. The CSV must contain at least the columns XXX.

The following parameters can be tuned for this step:

- XXX
- XXX
  
The available tools for this step are XXX.

**NB: When running this step, XXXX***

#### Examples

Minimal input file:

```bash
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_1.fastq.gz,test_2.fastq.gz
```

In this example, the sample comes from multiple patients:

```bash
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_L001_1.fastq.gz,test_L001_2.fastq.gz
patient1,test_sample,lane_2,test_L002_1.fastq.gz,test_L002_2.fastq.gz
patient1,test_sample,lane_3,test_L003_1.fastq.gz,test_L003_2.fastq.gz
```

### Clone Tree Inference

This step can be started either from XXX files or XXXX. The CSV must contain at least the columns XXX.

The following parameters can be tuned for this step:

- XXX
- XXX
  
The available tools for this step are XXX.

**NB: When running this step, XXXX***

#### Examples

Minimal input file:

```bash
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_1.fastq.gz,test_2.fastq.gz
```

In this example, the sample comes from multiple patients:

```bash
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_L001_1.fastq.gz,test_L001_2.fastq.gz
patient1,test_sample,lane_2,test_L002_1.fastq.gz,test_L002_2.fastq.gz
patient1,test_sample,lane_3,test_L003_1.fastq.gz,test_L003_2.fastq.gz
```

### Starting with Signature Deconvolution

This step can be started from `rds` multisample CNAqc object. The CSV must contain at least the columns:

```bash
dataset,joint_table

```
When running this step, the mutation's information of the dataset of interest is extracted from multisample `mCNAqc` object and converted to `txt` format in order to import the constructed data file in R or Python. Mutation frequency count data is generated from point mutation data-frames. The optimal signature number and sparsity is determined by cross-validation (SparseSignatures), followed by discovering the signatures within the dataset.

#### Examples

Minimal input file:

```bash
dataset,joint_table
CLL,mcnaqc_miltisample.rds
``` -->

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/tumourevo
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/tumourevo releases page](https://github.com/nf-core/tumourevo/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
    - A generic configuration profile to be used with [Docker](http://docker.com/)
    - Pulls software from dockerhub: [`nfcore/tumourevo`](http://hub.docker.com/r/nfcore/tumourevo/)
- `singularity`
    - A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    - Pulls software from DockerHub: [`nfcore/tumourevo`](http://hub.docker.com/r/nfcore/tumourevo/)
- `conda`
    - Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
    - A generic configuration profile to be used with [Conda](https://conda.io/docs/)
    - Pulls most software from [Bioconda](https://bioconda.github.io/)
- `test`
    - A profile with a complete configuration for automated testing
    - Includes links to test data so needs no other parameters

<!-- TODO nf-core: Document required command line parameters -->

<!-- ### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--single_end`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--single_end --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run. -->

<!-- ## Reference genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

<!-- TODO nf-core: Update reference genome example according to what is needed -->

<!-- ```nextflow
params {
  genomes {
    'GRCh37' {
      fasta   = '<path to the genome fasta file>' // Used if no star index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
``` -->

<!-- TODO nf-core: Describe reference path flags -->

<!-- ### `--fasta`

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta '[path to Fasta reference]'
```

### `--igenomes_ignore`

Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`. -->

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

<!-- ## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't. -->

## Other command line parameters

<!-- TODO nf-core: Describe any other command line flags here -->

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
