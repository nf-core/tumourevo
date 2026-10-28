process JOIN_CNAQC {
    tag "$meta.id"
    label "process_low"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-cnaqc%3A1.1.5--r44hdfd78af_0':
        'quay.io/biocontainers/r-cnaqc:1.1.5--r44hdfd78af_0' }"

    input:
    tuple val(meta), path(rds_list, stageAs: 'input*.rds'), val(tumour_samples)

    output:
    tuple val(meta), path("*ALL.rds"), val(tumour_samples),  emit: rds_all,  optional: true
    tuple val(meta), path("*PASS.rds"), val(tumour_samples), emit: rds_pass, optional: true
    path "versions.yml",                                     emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "$meta.id"
    def qc_filter = args.qc_filter != null ? args.qc_filter : false
    def keep_original = args!="" && args.keep_original ? "$args.keep_original" : ""
    def qc_chr = args.qc_chr != null ? args.qc_chr : false


    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(CNAqc)

    samples = substr("$tumour_samples", 2, nchar("$tumour_samples")-1)
    samples = strsplit(samples, ", ")[[1]]

    result = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(file){
                readRDS(file)
            })
    names(result) = samples


    if (as.logical("$qc_chr") == TRUE) {
        result = lapply(result, function(x) {
            x = lapply(x, function(chr) {
                chr\$mutations = chr\$mutations %>% dplyr::rename(Indiv = sample) %>% dplyr::select(-additional_info)
                return(chr)
            })
            return(x)
        })

        chromosomes = lapply(result, names) %>% unlist %>% unique

        out_all = lapply(chromosomes, function(cc) {
            mcnaqc_list = lapply(result, function(df) {
                df[[cc]]
            })

            out_all_by_chr = multisample_init(mcnaqc_list,
                                                QC_filter = FALSE,
                                                keep_original = as.logical("$keep_original"),
                                                discard_private = FALSE)
            return(out_all_by_chr)
        })
        names(out_all) = chromosomes

        saveRDS(object = out_all, file = paste0("$prefix", "_by_chr_multi_cnaqc_ALL.rds"))

        if (as.logical("$qc_filter") == TRUE){

          # select the chr elements
          out_PASS = lapply(chromosomes, function(cc) {
            mcnaqc_list = lapply(result, function(df) {
              df[[cc]]
            })
            # open the trycatch per chr
            tryCatch(expr = {
              out_PASS_by_chr = multisample_init(mcnaqc_list,
                                                        QC_filter = TRUE,
                                                        keep_original = as.logical("$keep_original"),
                                                        discard_private = FALSE)
              return(out_PASS_by_chr)
            }, error = function(e) {
              print(e)
              print('Not found common segments with QC PASS karyotype, the multi-CNaqc object will be NULL: re-run the pipeline with --filter false')
              out_PASS <<- NULL
            }
            )
          })
          names(out_PASS) = chromosomes
          saveRDS(object = out_PASS, file = paste0("$prefix", "_by_chr_multi_cnaqc_PASS.rds"))
        }

    } else {

        for (name in names(result)){
            result[[name]]\$mutations = result[[name]]\$mutations %>% dplyr::rename(Indiv = sample) %>% dplyr::select(-additional_info)
        }

        out_all = CNAqc::multisample_init(result,
                                QC_filter = FALSE,
                                keep_original = as.logical("$keep_original"),
                                discard_private = FALSE)

        saveRDS(object = out_all, file = paste0("$prefix", "_multi_cnaqc_ALL.rds"))


        if (as.logical("$qc_filter") == TRUE){
          tryCatch(expr = {
            out_PASS = CNAqc::multisample_init(result,
                                QC_filter = TRUE,
                                keep_original = as.logical("$keep_original"),
                                discard_private = FALSE)
            }, error = function(e) {
                print(e)
                print('Not found common segments with QC PASS karyotype, the multi-CNaqc object will be NULL: re-run the pipeline with --filter false')
                out_PASS <<- NULL
            }
          )
          saveRDS(object = out_PASS, file = paste0("$prefix", "_multi_cnaqc_PASS.rds"))
        }
    }

    # version export
    f <- file("versions.yml","w")
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    cnaqc_version <- sessionInfo()\$otherPkgs\$CNAqc\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    CNAqc:", cnaqc_version), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    close(f)

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_multi_cnaqc_ALL.rds
    touch ${prefix}_multi_cnaqc_PASS.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CNAqc: \$(Rscript -e "cat(as.character(packageVersion('CNAqc')))")
        dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
