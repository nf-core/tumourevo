process GET_POSITIONS_ALL {
    tag "$meta.id"
    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/913112a2d9295d35fe97874caf5f500df5a98d6fef1cb5861fd64caa0223a047/data':
        'community.wave.seqera.io/library/r-cnaqc_r-cli_r-dplyr_r-readr_pruned:0fc82bfd06afe6dc' }"

    input:
    tuple val(meta), path(rds_list, stageAs: "*.rds")

    output:
    tuple val(meta), path("*_all_positions.rds"), emit: all_pos
    path "versions.yml",                          emit: versions

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)

    positions = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(rds){
        df = readRDS(rds)

        df1 = df[[1]]\$mutations
        df2 = df[[2]]\$mutations
        if ('HGVSp' %in% colnames(df1)){
            df = df1
        } else {
            df = df2
        }

        df = df %>%
              dplyr::mutate(IMPACT = factor(IMPACT, levels = c("LOW", "MODIFIER", "MODERATE", "HIGH"))) %>%
              dplyr::mutate(pos = paste(chr, from, to, sep = ":")) %>%
              dplyr::group_by(pos) %>%
      	      dplyr::arrange(desc(SYMBOL)) %>%
              dplyr::slice_max(order_by = as.numeric(IMPACT), n = 1, with_ties = FALSE) %>%
	            dplyr::ungroup() %>%
              dplyr::select(chr, from, to, ref, alt, SYMBOL, HGVSp, IMPACT)
    })
    all = positions %>% dplyr::bind_rows() %>% dplyr::distinct()
    saveRDS(object = all, file = paste0("$prefix", "_all_positions.rds"))

    # version export
    f <- file("versions.yml","w")
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    close(f)

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_all_positions.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
