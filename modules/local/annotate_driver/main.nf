process ANNOTATE_DRIVER {
    tag "$meta.id"
    label "process_single"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/cnaqc:version1.0' :
        'docker.io/lvaleriani/cnaqc:version1.0' }"

    input:
    tuple val(meta), path(rds), path(driver_list)

    output:
    tuple val(meta), path("*.rds"),     emit: rds
    path "versions.yml",                emit: versions

    script:
    def args    =   task.ext.args   ?: ''
    def prefix  =   task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)
    library(tidyr)

    data = readRDS("$rds")
    SNV = data[["$meta.tumour_sample"]]
    SNV = SNV\$mutations

    drivers_table = readr::read_tsv(file = "$driver_list")

    if("$meta.cancer_type" %in% drivers_table\$CANCER_TYPE){
        drivers_table = drivers_table %>%
            dplyr::group_by(SYMBOL) %>%
            dplyr::reframe(CGC_CANCER_GENE = any(CGC_CANCER_GENE), dplyr::across(dplyr::everything())) %>%
            dplyr::filter(CGC_CANCER_GENE) %>%
            dplyr::filter(CANCER_TYPE == "$meta.cancer_type")
    } else {
        drivers_table = drivers_table %>%
            dplyr::group_by(SYMBOL) %>%
            dplyr::reframe(CGC_CANCER_GENE = any(CGC_CANCER_GENE), dplyr::across(dplyr::everything())) %>%
            dplyr::filter(CGC_CANCER_GENE) %>%
            dplyr::mutate(CANCER_TYPE = "PANCANCER")
    }

    drivers_table = drivers_table %>%
        dplyr::select(SYMBOL, CANCER_TYPE, CGC_CANCER_GENE) %>%
        unique()

    cancer_type = "$meta.cancer_type"
    if (!(cancer_type) %in% drivers_table\$CANCER_TYPE ){
        cancer_type = 'PANCANCER'
    }

    x = SNV %>%
        dplyr::mutate(CANCER_TYPE = cancer_type) %>%
        dplyr::left_join(
            drivers_table,
            by = c('SYMBOL', 'CANCER_TYPE')
        ) %>%
        tidyr::separate(HGVSp, ':', into = c('s1', 's2'), remove=F) %>%
        dplyr::mutate(tmp_s2 = ifelse(is.na(s2), '', paste0('_', s2))) %>%
        dplyr::mutate(
            is_driver = (CGC_CANCER_GENE & IMPACT %in% c('MODERATE', 'HIGH')),
            driver_label = paste0(SYMBOL, tmp_s2)
        ) %>%
        select(-tmp_s2) %>%
        mutate(is_driver = ifelse(is.na(is_driver), FALSE, is_driver))

    data[["$meta.tumour_sample"]]\$mutations = x
    saveRDS(object = data, file = paste0("$prefix", "_driver.rds"))

    # version export
    f <- file("versions.yml","w")
    readr_version <- sessionInfo()\$otherPkgs\$readr\$Version
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    tidyr_version <- sessionInfo()\$otherPkgs\$tidyr\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    readr:", readr_version), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    writeLines(paste("    tidyr:", tidyr_version), f)
    close(f)

    """
}
