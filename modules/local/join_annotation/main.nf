
process JOIN_ANNOTATION {
    tag "$meta.id"
    container = 'docker://lvaleriani/cnaqc:dev1'

    input:

    tuple val(meta), path(annot_rds), val(meta2), path(dnds_rds) 

    output:
    tuple val(meta), path("*.rds"), emit: rds


    script:

    def args = task.ext.args ?: ''    
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)

    data = readRDS("$annot_rds")
    SNV = data[["$meta.tumour_sample"]] 
    annot = SNV\$mutations

    dnds = readRDS("$dnds_rds")
    dnds_driver = dnds %>% select(chr, from, potential_driver)

    annot = left_join(annot, dnds_driver) %>% mutate(is_driver = ifelse(know_driver == TRUE | potential_driver == TRUE, TRUE, FALSE))
    data[["$meta.tumour_sample"]]\$mutations = annot
    saveRDS(data, paste0("$prefix", "_driver_annotation.rds"))
    
    """
}
