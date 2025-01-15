process GET_POSITIONS_REL {
    tag "$meta.id"
    container = 'docker://lvaleriani/cnaqc:dev1'

    input:
    tuple val(meta), path(rds), path(all_pos)

    output:
    tuple val(meta), path("*_positions_missing"),   emit: bed
    path "versions.yml",                            emit: versions


    script:
    def args    = task.ext.args     ?:  ''
    def prefix  = task.ext.prefix   ?:  "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)

    all_positions = readRDS("$all_pos") %>%
                dplyr::mutate(id = paste(chr, from, to, sep = ":"))

    df = readRDS("$rds")
    positions = df[["$meta.tumour_sample"]]\$mutations %>%
                dplyr::mutate(id = paste(chr, from, to, sep = ":")) %>%
                dplyr::select(chr, from, to, ref, alt, id)

    missed = all_positions %>% filter(!(id %in%  positions\$id)) %>%
            dplyr::filter(chr %in% c(paste0('chr', seq(1,22)), 'chrX', 'chrY')) %>%
            dplyr::select(chr, from)

    write.table(file = paste0("$meta.tumour_sample", "_positions_missing"), missed, quote = F, sep = "\t", row.names = F, col.names = F)

    # version export
    f <- file("versions.yml","w")
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    close(f)

    """
}
