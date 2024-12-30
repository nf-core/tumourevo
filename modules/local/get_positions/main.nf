process GET_POSITIONS_ALL {
    tag "$meta.id"
    container = 'docker://lvaleriani/cnaqc:dev1'

    input:
    tuple val(meta), path(rds_list, stageAs: '*.rds')

    output:

    tuple val(meta), path("*_all_positions.rds"), emit: all_pos

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    positions = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(rds){
        df = readRDS(rds)
        df = df[[1]]\$mutations %>% dplyr::mutate(id = paste(chr, from, to, sep = ":")) %>% dplyr::select(chr, from, to, ref, alt, id)
    })
    all = positions %>% dplyr::bind_rows() %>% dplyr::distinct() %>% dplyr::select(-id)
    saveRDS(object = all, file = paste0("$prefix", "_all_positions.rds"))
    #write.table(file = paste0("$prefix", "_all_positions.bed"), all, quote = F, sep = "\t", row.names = F, col.names = T)
    """
}
