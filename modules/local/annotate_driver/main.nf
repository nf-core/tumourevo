process ANNOTATE_DRIVER {
    tag "$meta.id"
    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/913112a2d9295d35fe97874caf5f500df5a98d6fef1cb5861fd64caa0223a047/data':
        'community.wave.seqera.io/library/r-cnaqc_r-cli_r-dplyr_r-readr_pruned:0fc82bfd06afe6dc' }"

    input:
    tuple val(meta), path(rds), path(driver_list)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    path "versions.yml",            emit: versions

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)
    library(tidyr)

    data = readRDS("$rds")
    SNV = data[["$meta.tumour_sample"]]
    SNV = SNV\$mutations

    drivers_table = readr::read_tsv(file = "$driver_list")

    tumour_type = "$meta.cancer_type"
    if(tumour_type %in% drivers_table\$TUMOUR_TYPE){
        drivers_table = drivers_table %>%
            dplyr::filter(TUMOUR_TYPE == tumour_type)
    } else {
        drivers_table = drivers_table %>%
            dplyr::mutate(TUMOUR_TYPE = "PANCANCER")
        tumour_type = 'PANCANCER'
    }

    drivers_table = drivers_table %>%
        dplyr::select(SYMBOL, TUMOUR_TYPE) %>%
        dplyr::distinct()

    x = SNV %>%
        dplyr::left_join(
            drivers_table,
            by = c('SYMBOL')
        ) %>%
        tidyr::separate(HGVSp, ':', into = c('s1', 's2'), remove=F) %>%
        dplyr::mutate(tmp_s2 = ifelse(is.na(s2), '', paste0('_', s2))) %>%
        dplyr::mutate(
            is_driver = (TUMOUR_TYPE != "" & SYMBOL != "" & IMPACT %in% c('MODERATE', 'HIGH')),
            driver_label = paste0(SYMBOL, tmp_s2)
        ) %>%
        dplyr::select(-tmp_s2) %>%
        dplyr::mutate(is_driver = ifelse(is.na(is_driver), FALSE, is_driver)) %>%
        dplyr::mutate(IMPACT = factor(IMPACT, levels = c("LOW", "MODIFIER", "MODERATE", "HIGH"))) %>%
        dplyr::mutate(pos = paste(chr, from, to, sep = ":")) %>%
        dplyr::group_by(pos) %>%
        dplyr::mutate(additional_info = tibble(SYMBOL, IMPACT, Consequence, HGVSc, HGVSp, is_driver)) %>%
        dplyr::mutate(additional_info = list(across(additional_info))) %>%
	      dplyr::arrange(desc(SYMBOL)) %>%
        dplyr::slice_max(order_by = as.numeric(IMPACT), n = 1, with_ties = FALSE) %>%
        dplyr::distinct(pos, SYMBOL, IMPACT, .keep_all = TRUE) %>%
        dplyr::ungroup()

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

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_driver.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        tidyr: \$(Rscript -e "cat(as.character(packageVersion('tidyr')))")
    END_VERSIONS
    """
}
