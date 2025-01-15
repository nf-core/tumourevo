process VCF_PROCESSING {
    tag "$meta.id"

    container "docker://lvaleriani/cnaqc:version1.0"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.rds"),     emit: rds
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                = task.ext.args         ?: ''
    def prefix                  = task.ext.prefix   ?: "${meta.id}"
    def filter_mutations    = args.filter_mutations ? "$args.filter_mutations": "NULL"

    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(tidyr)
    library(vcfR)

    source(paste0("$moduleDir", '/parser_vcf.R'))

    # Read vcf file
    vcf = vcfR::read.vcfR("$vcf")

    # Check from which caller the .vcf has been produced
    source = vcfR::queryMETA(vcf, element = 'source')[[1]]

    if (TRUE %in% grepl(pattern = 'Mutect', x = source)){
        calls = parse_Mutect(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample", filter_mutations = as.logical("$filter_mutations"))

    } else if (TRUE %in% grepl(pattern = 'strelka', x = source)){
        calls = parse_Strelka(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample", filter_mutations = as.logical("$filter_mutations"))

    } else if (TRUE %in% grepl(pattern = 'Platypus', x = source)){
        calls = parse_Platypus(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample", filter_mutations = as.logical("$filter_mutations"))

    } else if (TRUE %in% grepl(pattern = 'freeBayes', x = source)){
        calls = parse_FreeBayes(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample", filter_mutations = as.logical("$filter_mutations"))

    } else {
        stop('Variant Caller not supported.')
    }

    saveRDS(object = calls, file = paste0("$prefix", "_snv.rds"))

    # version export
    f <- file("versions.yml","w")
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    tidyr_version <- sessionInfo()\$otherPkgs\$tidyr\$Version
    vcfR_version <- sessionInfo()\$otherPkgs\$vcfR\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    writeLines(paste("    tidyr:", tidyr_version), f)
    writeLines(paste("    vcfR:", vcfR_version), f)
    close(f)

    """
}

