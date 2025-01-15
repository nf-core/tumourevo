process CNA_PROCESSING {
    tag "$meta.id"

    container "docker://lvaleriani/cnaqc:version1.0"

    input:
    tuple val(meta), path(cna_segs), path(cna_extra)

    output:
    tuple val(meta), path("*.rds"),                            emit: rds
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?:  ''
    def prefix  = task.ext.prefix   ?:  "${meta.id}"

    """
    #!/usr/bin/env Rscript
    source(paste0("$moduleDir", '/parser_CNA.R'))

    if ("$meta.cna_caller" == 'sequenza'){
    CNA = parse_Sequenza(segments = "$cna_segs", extra = "$cna_extra")

    } else if ("$meta.cna_caller" == 'ASCAT'){
    CNA = parse_ASCAT(segments = "$cna_segs", extra = "$cna_extra")

    } else {
    stop('Copy Number Caller not supported.')
    }

    saveRDS(object = CNA, file = paste0("$prefix", "_cna.rds"))

    # version export
    f <- file("versions.yml","w")
    readr_version <- sessionInfo()\$otherPkgs\$readr\$Version
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    readr:", readr_version), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    close(f)

    """
}
