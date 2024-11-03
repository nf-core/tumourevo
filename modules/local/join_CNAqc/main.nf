process JOIN_CNAQC {
  tag "$meta.id"
  container = 'docker://lvaleriani/cnaqc:dev1'

  input:
    
    tuple val(meta), path(rds_list, stageAs: '*.rds'), val(tumour_samples)
  
  output:

    tuple val(meta), path("*.rds"), val(tumour_samples), emit: rds, optional: true

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
      
    out = CNAqc::multisample_init(result, 
                              QC_filter = as.logical("$qc_filter"), 
                              keep_original = as.logical("$keep_original"), 
                              discard_private = FALSE)
    saveRDS(object = out, file = paste0("$prefix", "_multi_cnaqc.rds"))
    """
}
