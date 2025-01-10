//
// Mutations extraction from mCNAqc
//

process RDS_PROCESSING {
    tag "$meta.id"
    container = 'docker://lvaleriani/cnaqc:version1.0'

    input:
        tuple val(meta), path(rds_join), val(tumour_samples)

    output:
        tuple val(meta), path("*_joint_table.tsv"), val(tumour_samples), emit: tsv

    script:
    def args                              = task.ext.args                                 ?: ''
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(CNAqc)
    library(tidyr)

    source("$moduleDir/utils.R")


    multi_cnaqc = readRDS(file = "$rds_join")
    mutations_multisample <- get_sample(m_cnaqc_obj = multi_cnaqc,sample = get_sample_name(multi_cnaqc),
                                    which_obj = "original")
    multisample_jointTable <- list()

    for (s in get_sample_name(multi_cnaqc)){
        purity <- mutations_multisample[[s]][["purity"]]
        multisample_jointTable[[s]] <- mutations_multisample[[s]][["mutations"]] %>%
            dplyr::mutate(purity=purity)
    }

    joint_table <- bind_rows(multisample_jointTable)

    write.table(joint_table, file = paste0("$prefix","_joint_table.tsv"), append = F, quote = F, sep = "\t", row.names = FALSE)

    """
}
