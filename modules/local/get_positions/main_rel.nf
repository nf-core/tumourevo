process GET_POSITIONS_REL {
    tag "$meta.id"
    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/913112a2d9295d35fe97874caf5f500df5a98d6fef1cb5861fd64caa0223a047/data':
        'community.wave.seqera.io/library/r-cnaqc_r-cli_r-dplyr_r-readr_pruned:0fc82bfd06afe6dc' }"

    input:
    tuple val(meta), path(rds), path(all_pos)

    output:
    tuple val(meta), path("*_positions_missing.bed"), emit: bed
    path "versions.yml",                          emit: versions

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)

    all_positions = readRDS("$all_pos") %>%
                dplyr::mutate(id = paste(chr, from, to, sep = ":"))

    df = readRDS("$rds")
    positions = df[["$meta.tumour_sample"]]\$mutations %>%
                dplyr::mutate(id = paste(chr, from, to, sep = ":")) %>%
                dplyr::select(chr, from, to, ref, alt, id)

    missed = all_positions %>%
            dplyr::filter(!(id %in%  positions\$id)) %>%
            dplyr::filter(chr %in% c(paste0('chr', seq(1,22)), 'chrX', 'chrY')) %>%
            dplyr::select(chr, from, to) %>%
            dplyr::mutate(from = from-1, to = to-1)

    write.table(file = paste0("$meta.id", "_positions_missing.bed"), missed, quote = F, sep = "\t", row.names = F, col.names = F)

    # version export
    f <- file("versions.yml","w")
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    close(f)

    """

    stub:
    """
    touch ${meta.tumour_sample}_positions_missing.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
