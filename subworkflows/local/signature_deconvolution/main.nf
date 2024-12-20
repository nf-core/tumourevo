//
// MUTATIONAL SIGNATURES DECONVOLUTION WORKFLOW
//

include { FORMATTER as FORMATTER_RDS_SIGPROFILER } from "../../../subworkflows/local/formatter/main"
include { FORMATTER as FORMATTER_RDS_SPARSESIGNATURES } from "../../../subworkflows/local/formatter/main"
include { SPARSE_SIGNATURES } from "../../../modules/local/SparseSignatures/main"
include { DOWNLOAD_GENOME_SIGPROFILER } from "../../../modules/local/SigProfiler/download/main"
include { SIGPROFILER } from "../../../modules/local/SigProfiler/SigProfiler/main"


workflow SIGNATURE_DECONVOLUTION {
    take: 
    join_cnaqc_out // tuple val(meta), path("*.rds"), val(tumor_samples)

    main:
    plot_pdf = null
    plot_rds = null
    signatures_nmfOut = null 
    bestConf = null
    sign_cv = null
    genome_path = null
    Sigprofiler_out = null


    if (params.tools && params.tools.split(',').contains('sparsesignatures')) {
        FORMATTER_RDS_SPARSESIGNATURES(join_cnaqc_out, "rds")
        input_sparsesig = FORMATTER_RDS_SPARSESIGNATURES.out.map { meta, tsv, sample ->
            meta = meta + [id: "${meta.dataset}"] 
            [meta.subMap('dataset', 'id'), tsv] }
            | groupTuple
  
        SPARSE_SIGNATURES(input_sparsesig) // run SparseSignatures
        
        plot_pdf = SPARSE_SIGNATURES.out.signatures_plot_pdf
        plot_rds = SPARSE_SIGNATURES.out.signatures_plot_rds
        signatures_nmfOut = SPARSE_SIGNATURES.out.signatures_nmfOut_rds
        bestConf = SPARSE_SIGNATURES.out.signatures_bestConf_rds
        sign_cv = SPARSE_SIGNATURES.out.signatures_cv_rds
    } 
  

    if (params.tools && params.tools.split(',').contains('sigprofiler')) {
      
        if (params.download_sigprofiler_genome) {    
            genome_path = DOWNLOAD_GENOME_SIGPROFILER(params.genome).genome_sigprofiler           
        } else {
       
            genome_path = params.genome_installed_path
        }
        
        FORMATTER_RDS_SIGPROFILER(join_cnaqc_out, "rds")

        input_sigprofiler = FORMATTER_RDS_SIGPROFILER.out.map { meta, tsv, sample ->
            meta = meta + [id: "${meta.dataset}"] 
            [meta.subMap('dataset', 'id'), tsv]}
            | groupTuple  
  
        Sigprofiler_out = SIGPROFILER(input_sigprofiler, genome_path)
    }

    emit:
    plot_pdf
    plot_rds
    signatures_nmfOut
    bestConf
    sign_cv
    Sigprofiler_out
    
}
