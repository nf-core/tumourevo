//
// DNDSCV PROCESS
//
process DNDSCV {
  debug true
  tag "$meta.id"
  container "${workflow.containerEngine == 'singularity' ? 'docker://tucano/dndscv:latest' : 'tucano/dndscv:latest'}"


  input:
    
    tuple val(meta), path(rds), val(patients), val(samples), path(reference)
  
  output:

    tuple val(meta), path("*_annot.rds"), emit: annot_rds
    tuple val(meta), path("*_dnds.rds"), emit: dnds_rds
  
  script:

    def args    = task.ext.args ?: ""
    def prefix  = task.ext.prefix ?: "${meta.id}"  
    
    """
    #!/usr/bin/env Rscript
    library(stringr)
    library(readr)
    library(dndscv)
    library(dplyr)


    # read list of rds
    rds_list <- strsplit("$rds", " ")[[1]]
    samples <- strsplit(stringr::str_sub("$samples", 2, -2),", ")[[1]]

    mut_list <- list()
    for (idx in 1:length(samples)){
      rds <- readRDS(rds_list[idx])
      mut <- rds[[samples[idx]]]\$mutations %>% mutate(sample_id = rep(samples[idx], n()))
      mut_list[[samples[idx]]] <- mut
    }
    mutations <- bind_rows(mut_list) %>% 
      select(chr, from, to, ref, alt, sample_id)

    reference <- "$reference"
    load(reference)

    chr_ref <- as.character(RefCDS[[1]][['chr']])
    reference_with_chr <- startsWith(chr_ref, "chr")
    chr_mut <- as.character(mutations[['chr']][1])
    mutations_with_chr <- startsWith(chr_mut,"chr")

    dndscv_input <- mutations |>
      dplyr::select(chr,from,ref,alt) |>
      unique() |>
      dplyr::rename(pos="from") |>
      dplyr::mutate(sample="$prefix", .before = chr) |>
      dplyr::mutate(chr=str_replace(chr,"chr",""))

    # add chr in dndscv_input when needed
    if (reference_with_chr) {
      dndscv_input <- dndscv_input |> dplyr::mutate(chr=paste("chr",chr,sep=""))
    }

    dndscv_result <- dndscv::dndscv(
      mutations = dndscv_input,
      outmats = TRUE,
      max_muts_per_gene_per_sample = Inf,
      max_coding_muts_per_sample = Inf,
      outp = 3,
      use_indel_sites = TRUE,
      min_indels = 1,
      refdb=reference,
      cv=NULL
    )
    
    sel <- left_join(dndscv_result[['sel_cv']], dndscv_result[['sel_loc']])
    annotation <- left_join(dndscv_result\$annotmuts, sel, by = c("gene"="gene_name"))
    annotation <- annotation |> 
      dplyr::select(!sampleID) |> 
      dplyr::rename(from="pos", alt="mut")
    
    annotation <- annotation |> 
      dplyr::mutate(
        potential_driver = (qmis_cv <= 0.1 | qtrunc_cv <= 0.1 | qallsubs_cv <= 0.1)
      )

    mutations <- mutations |> dplyr::mutate(chr=str_replace(chr,"chr",""))
    annotation <- annotation |> dplyr::mutate(chr=str_replace(chr,"chr",""))
    
    if (mutations_with_chr) {
      mutations <- mutations |> dplyr::mutate(chr=paste("chr",chr,sep=""))
      annotation <- annotation |> dplyr::mutate(chr=paste("chr",chr,sep=""))
    }

    saveRDS(annotation, paste0("$prefix", "_annot.rds"))
    saveRDS(dndscv_result, paste0("$prefix", "_dnds.rds"))

  """
}


// dndscv_runner.R \\
//   -i ${rds} \\
//   -s ${samples} \\
//   -r ${reference} \\
//   -o "${prefix}_dnds.rds" \\
//   --qmis_cv ${params.dndscv_qmis_cv} \\
//   --qtrunc_cv ${params.dndscv_qtrunc_cv} \\
//   --qallsubs_cv ${params.dndscv_qallsubs_cv} \\
//   ${args}
