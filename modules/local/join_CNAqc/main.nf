process JOIN_CNAQC {
  tag "$meta.id"
  container = 'docker://lvaleriani/cnaqc:dev1'

  input:
    
    tuple val(meta), path(rds_list, stageAs: '*.rds'), val(tumour_samples)
  
  output:

    tuple val(meta), path("*ALL.rds"), val(tumour_samples), emit: rds_all, optional: true
    tuple val(meta), path("*PASS.rds"), val(tumour_samples), emit: rds_pass, optional: true

  script:

    def args                                = task.ext.args    
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"                                      ?: ''
    def qc_filter                           = args!='' && args.qc_filter                            ?  "$args.qc_filter" : ""
    def keep_original                       = args!='' && args.keep_original                        ?  "$args.keep_original" : ""

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(CNAqc)

    samples = substr("$tumour_samples", 2, nchar("$tumour_samples")-1)
    samples = strsplit(samples, ", ")[[1]]

    result = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(file){
              readRDS(file)
              }) 
    names(result) = samples
      
    for (name in names(result)){
        result[[name]]\$mutations = result[[name]]\$mutations %>% dplyr::rename(Indiv = sample)
    }
      
    out_all = CNAqc::multisample_init(result, 
                              QC_filter = FALSE, 
                              keep_original = as.logical("$keep_original"), 
                              discard_private = FALSE)

    out_PASS = CNAqc::multisample_init(result, 
                              QC_filter = TRUE, 
                              keep_original = as.logical("$keep_original"), 
                              discard_private = FALSE)

    saveRDS(object = out_all, file = paste0("$prefix", "_multi_cnaqc_ALL.rds"))
    saveRDS(object = out_PASS, file = paste0("$prefix", "_multi_cnaqc_PASS.rds"))
    """
}
