process GET_POSITIONS_ALL {
    tag "$meta.id"
    container = 'docker://lvaleriani/cnaqc:version1.0'

    input:
    tuple val(meta), path(rds_list, stageAs: '*.rds')

    output:
    tuple val(meta), path("*_all_positions.rds"),   emit: all_pos
    path "versions.yml",                            emit: versions


    script:
    def args    = task.ext.args     ?:  ''
    def prefix  = task.ext.prefix   ?:  "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)

    positions = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(rds){
        df = readRDS(rds)
        df = df[[1]]\$mutations %>% dplyr::mutate(id = paste(chr, from, to, sep = ":")) %>% dplyr::select(chr, from, to, ref, alt, id)
    })
    all = positions %>% dplyr::bind_rows() %>% dplyr::distinct() %>% dplyr::select(-id)
    saveRDS(object = all, file = paste0("$prefix", "_all_positions.rds"))

    # version export
    f <- file("versions.yml","w")
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    close(f)

    """
}
