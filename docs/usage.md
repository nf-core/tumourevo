# nf-core/tumourevo: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/tumourevo/usage](https://nf-co.re/tumourevo/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**tumourevo** is a workflow to infer a tumour evolution model from whole-genome sequencing (WGS) data.

Through the analysis of variant and copy-number calls, it reconstructs the evolutionary process leading to the observed tumour genome. Most of the analyses can be done at mutliple levels: single sample, multiple samples from the same patient (multi-region/longitudinal assays), and multiple patients from distinct cohorts.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

/_The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet.
The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to
match those defined in the table below.
_/

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the parameter `--input` to specify its location. It has to be a comma-separated file with at least 5 columns, and a header row as shown in the examples below.

It is recommended to use the absolute path of the files, but a relative path should also work.

For the joint analysis of multiple samples, a tumour BAM file is required for each sample, such that the number of reads of a private mutation can be retrieved for all the samples thorugh `mpileup`.

Multiple samples from the same patient must be specified with the same `dataset` ID, `patient` ID, and a different `tumour_sample` ID. `normal_sample` columns is required.

Multiple patients from the same dataset must be specified with the same `dataset` ID, and a different `patient` ID.

**tumourevo** will output sample-specific results in a different directory for _each sample_, patient-specific results in a common directory for _each patient_, and dataset-specific results in a common directory for _each dataset_.

Output from different workflows, subworkflows and modules will be in a specific directory for each dataset, patient, sample and tool configuration.

A minimal input sample sheet example for two samples from the same patient:

```csv title="samplesheet.csv"
dataset,patient,sample,normal_sample,vcf,vcf_tbi,cna_segments,cna_extra,cna_caller,cancer_type
dataset1,patient1,S1,N1,patient1_S1.vcf.gz,patient1_S1.vcf.gz.tbi,/CNA/patient1/S1/segments.txt,/CNA/patient1/S1/purity_ploidy.txt,caller,PANCANCER
dataset1,patient1,S2,N1,patient1_S1.vcf.gz,patient1_S1.vcf.gz.tbi,/CNA/patient1/S2/segments.txt,/CNA/patient1/S2/purity_ploidy.txt,caller,PANCANCER
```

#### Overview: Samplesheet Columns

| Column                   | Description                                                                                                                                                                                                                               |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --- |
| `dataset`                | **Dataset ID**; when sequencing data from multiple datasets is analysed, it designates the source dataset of each patient; must be unique for each dataset, but one dataset can contain samples from multiple patients. <br /> _Required_ |
| `patient`                | **Patient ID**; designates the patient/subject; must be unique for each patient, but one patient can have multiple samples (e.g. from multiple regions or multiple time points). <br /> _Required_                                        |
| `tumour_sample`          | **Sample ID** for each sample; more than one sample for each subject is possible. Must match the sample ID present in the VCF. <br /> _Required_                                                                                          |
| `normal_sample`          | **Normal sample ID** of each sample. Must match the normal sample ID present in the VCF. <br /> _Required_                                                                                                                                |     |
| `vcf`                    | Full path to the vcf file. <br /> _Required_                                                                                                                                                                                              |
| `tbi`                    | Full path to the vcf `tabix` index file. <br /> _Required_                                                                                                                                                                                |
| `cna_caller`             | Name of the copy number caller used to generate your data. <br /> _Required_                                                                                                                                                              |
| `cna_segments`           | Full path to the segmentation files and copy number state from copy-number calling. <br /> _Required_                                                                                                                                     |
| `cna_extra`              | Full path to files including the ploidy and purity estimate from the copy-number caller. <br /> _Required_                                                                                                                                |
| `cancer_type`            | Tumour type (either `PANCANCER` or one of the tumour type present in the driver table) <br /> _Required_                                                                                                                                  |
| `tumour_alignment`       | Full path to the tumour bam file. <br /> _Optional_                                                                                                                                                                                       |
| `tumour_alignment_index` | Full path to the tumour bam index file. <br /> _Optional_                                                                                                                                                                                 |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Pipeline modalities

The tumourevo pipeline supports variant annotation, driver annotation, quality control processes, subclonal deconvolution and signature deconvolution analysis through various tools. It can be used to analyse both single sample experiments and longitudinal/multi-region assays, in which multiple samples of the same patient are avaiable.
As input, you must provide at least information on the samples, the VCF file from one of the supported callers and the output of one of the supported copy number caller. By default, if multiple samples from the same patient are provided, they will be analysed in a multivariate framework (which affects in particular the subclonal deconvolution deconvolution steps) to retrieve information useful in the reconstruction of the evolutionary process. Depending on the variant calling strategy (single sample or multi sample) and the provided input files, different strategies will be applied.

<!-- aggiungi un riassunto di cosa voglia dire single e multi sample (analisi multivariata, soprattutto per subclonal deconv)
E' possibile usarla sia nel caso di vc multi sample che indipendente -->

#### Variant calling

##### 1. Multi-sample variant calling

Modern tools (ie: Platypus and Mutect2) allow to perform variant calling directly in multisample mode. If the VCFs provided as input are already multisample, no additional step is required.

<!-- If the variant calling had been performed indepentently on each sample from the same patient,

If you run the pipeline in `singlesample` mode, all the samples, even if belonging to the same patient, are assumed to be independent. In this framework, the subclonal deconvolution is affected, identifying clonal and subclonal composition of the sample starting from allele frequency of detected somatic variants and identifying mutagenic processes for each independent element.  -->

<!-- rephrase better -->
<!-- Si assume che i campioni siano indipendenti e che quindi la subclonal deconv viene svolta cercando popolazioni sottoclonali a lv di singolo campione. Sign deconv si cercano i processi mutagenici comuni in un dataset fatto di elementi indipenti -->

###### Examples

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

##### 2. Single sample variant calling

If the variant calling performed independently on each sample, even if coming from the same patient, you can provide the BAM and BAI files from each tumour sample. In this way, a classical pileup strategy will be used in order to retrieve the depth for all samples of private mutations, in order to correctly perform the subclonal deconvolution analysis.

Input file for two patients without joint variant calling, bam files available:

```bash
dataset,patient,sample,normal_sample,vcf,vcf_tbi,cna_segments,cna_extra,cna_caller,cancer_type,tumour_bam,tumour_bai
dataset1,patient1,S1,N1,patient1_S1.vcf.gz,patient1_S1.vcf.gz.tbi,/CNA/patient1/S1/segments.txt,/CNA/patient1/S1/purity_ploidy.txt,caller,PANCANCER,patient1/BAM/S1.bam,patient1/BAM/S1.bam.bai
dataset1,patient1,S2,N1,patient1_S2.vcf.gz,patient1_S2.vcf.gz.tbi,/CNA/patient1/S2/segments.txt,/CNA/patient1/S2/purity_ploidy.txt,caller,PANCANCER,,patient1/BAM/S2.bam,patient1/BAM/S2.bam.bai
```

If you can not include the bam files in the input csv, the pipeline will run anyway, treating each sample as independent.

<!-- You can use the `multisample` mode of tumourevo to analyse samples from multi-region and longitudinal assays. This allows you to track in space and time the existing tumour populations, and better understand its heterogeneity. This modality integrates data across multiple samples, thus improving the resolution of subclonal structures and providing insights into the evolutionary dynamics and progression of the tumour.
Two of the avaiable tools for subclonal deconvolution, `pyclonevi` and `viber` can by-design be run in multi-sample mode, inferring the subclonal structure of samples. If you add `mobster` to the `--tool` parameter when running the pipeline in this modality, it will be run at first on each individual sample (since the tool does not support at the moment multi-sample analysis) in order to recognize neutral tail mutations and remove them. The mutations data manipulated in this way will then be processed by either `pyclone`, `viber` or both using the multivariate subclonal deconvolution as described before.  -->

#### 3. Filtering data

During the QC step, the pipeline will combine purity, copy number and mutation data to perform quality control on the copy number calls, by applying the [CNAqc algorithm](https://caravagnalab.github.io/CNAqc/). Each segment (for each sample) will be flagged as passing or not the QC, in the given combination of estimated purity and ploidy. According to CNAqc, a badly called segment should be recalled with a different purity estimation, in order to obtain more reliable results.
After CNAqc quality control, all the segments (coming from samples of the same patient) are used to build a multi-sample CNAqc object, in which a common segmentation is applied. In this way, only those regions that are shared among all samples will be kept in the new object. It is possible to control wheter to include or not in the new segmentation the segments, for each sample, that do not pass the QC test using the `--filter` flag. If it is setted to true, only QC passing segments for each sample will be used to build the mCNAqc object and will then be passed to the subclonal deconvolution steps. This will lead to exclude some regions of the genome and the mutations that sit on it, but should result in more precise analyses. Otherwise, keeping also the segments that do not pass the QC will result in not loosing any mutation but might lead to less precise results in the subclonal deconvolution steps.

#### 4. Driver annotation

You can retrieve tumour-specific drivers in the driver annotation step by specifying the tumour type in the input csv. Pan-cancer drivers will be retrieved by specifying `PANCANCER` as tumour type in the input csv file.
For this step, we currently refer to [IntOGen latest release](https://www.nature.com/articles/s41568-020-0290-x), but it is also possible to provide a custom driver table that will be used in the analysis.
Please note that the tumour types reported in the input file must correspond to those present in the table used for the annotation (default driver table used can be found [here](https://github.com/caravagnalab/nextflow_modules/blob/main/2023-05-31_IntOGen-Drivers/Unfiltered_drivers.tsv))

### 5. Available tools

We report the different tools included in the pipeline.

1. **Gene annotation**

   - [EnsemblVEP]()

2. **Driver annotation**

3. **Quality control**

   - [TINC](https://caravagnalab.github.io/TINC/)
   - [CNAqc](https://caravagnalab.github.io/CNAqc/)

  <!-- #### TINC

  #### INCOMMON -->

1. **Subclonal deconvolution**

   - [MOBSTER](https://caravagnalab.github.io/mobster/)
   - [PyClone-VI](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03919-2)
   - [VIBER](https://caravagnalab.github.io/VIBER/index.html)
   - [Ctree](https://caravagnalab.github.io/ctree/)

2. **Signature deconvolution**

   - [SparseSignatures](https://github.com/danro9685/SparseSignatures)
   - [SigProfiler](https://cancer.sanger.ac.uk/signatures/tools/)

## Running the pipeline

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

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command,
you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/tumourevo \
 -r <VERSION> \
 -profile <PROFILE> \
 -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh38'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/tumourevo
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/tumourevo releases page](https://github.com/nf-core/tumourevo/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is suported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
