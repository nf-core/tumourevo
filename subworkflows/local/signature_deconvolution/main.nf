//
// MUTATIONAL SIGNATURES DECONVOLUTION WORKFLOW
//

include { FORMATTER as FORMATTER_RDS_SIGPROFILER } from "../../../subworkflows/local/formatter/main"
include { FORMATTER as FORMATTER_RDS_SPARSESIGNATURES } from "../../../subworkflows/local/formatter/main"
include { SPARSE_SIGNATURES } from "../../../modules/nf-core/sparsesignatures/main"
include { SIGPROFILER } from "../../../modules/nf-core/sigprofiler/main"


workflow SIGNATURE_DECONVOLUTION {
    take:
    join_cnaqc_out // tuple val(meta), path("*.rds"), val(tumor_samples)

    main:
    ch_versions = Channel.empty()
    plot_pdf = null
    plot_rds = null
    signatures_nmfOut = null
    bestConf = null
    sign_cv = null
    mut_counts = null
    genome_path = null
    sigprofiler_out = null

    if (params.tools && params.tools.split(',').contains('sparsesignatures')) {
        FORMATTER_RDS_SPARSESIGNATURES(join_cnaqc_out, "rds")
        ch_versions = ch_versions.mix(FORMATTER_RDS_SPARSESIGNATURES.out.versions)

        input_sparsesig = FORMATTER_RDS_SPARSESIGNATURES.out.out_data.map { meta, tsv, sample ->
            meta = meta + [id: "${meta.dataset}"]
            [meta.subMap('dataset', 'id'), tsv] }
            .groupTuple()

        SPARSE_SIGNATURES(input_sparsesig, params.genome) // run SparseSignatures

        ch_versions = ch_versions.mix(SPARSE_SIGNATURES.out.versions_sparsesignatures)
        plot_pdf = SPARSE_SIGNATURES.out.signatures_plot_pdf
        plot_rds = SPARSE_SIGNATURES.out.signatures_plot_rds
        signatures_nmfOut = SPARSE_SIGNATURES.out.signatures_nmfOut_rds
        bestConf = SPARSE_SIGNATURES.out.signatures_bestConf_rds
        sign_cv = SPARSE_SIGNATURES.out.signatures_cv_rds
        mut_counts = SPARSE_SIGNATURES.out.signatures_mutCounts_rds
    }


    if (params.tools && params.tools.split(',').contains('sigprofiler')) {
        if (params.download_sigprofiler_genome) {
            genome_path = channel.fromPath('opt/null')
        } else {
            genome_path = params.genome_installed_path
        }

        FORMATTER_RDS_SIGPROFILER(join_cnaqc_out, "rds")
        ch_versions = ch_versions.mix(FORMATTER_RDS_SIGPROFILER.out.versions)

        input_sigprofiler = FORMATTER_RDS_SIGPROFILER.out.out_data.map { meta, tsv, sample ->
            meta = meta + [id: "${meta.dataset}"]
            [meta.subMap('dataset', 'id'), tsv]}
            .groupTuple()

        ch_seed = channel.of(
                    "\tSeed",
                    "0\t315836613344651455423314970566641622150"
                ).collectFile(name: 'seed.txt', newLine: true, sort: false)

        SIGPROFILER(input_sigprofiler, params.genome, genome_path, ch_seed)
        ch_versions = ch_versions.mix(SIGPROFILER.out.versions)
        sigprofiler_out = SIGPROFILER.out.results_sigprofiler
    }

    emit:
    plot_pdf
    plot_rds
    signatures_nmfOut
    bestConf
    sign_cv
    mut_counts
    sigprofiler_out
    versions = ch_versions

}
