process SPARSE_SIGNATURES {
    tag "$meta.id"
    label "process_low_long"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/sparsesignature:version1.0' :
        'docker.io/lvaleriani/sparsesignature:version1.0' }"

    input:
        tuple val(meta), path(tsv_join,  stageAs: '*.tsv')

    output:
        tuple val(meta), path("*_cv_means_mse.rds"), emit: signatures_cv_rds
        tuple val(meta), path("*_best_params_config.rds"), emit: signatures_bestConf_rds
        tuple val(meta), path("*_nmf_Lasso_out.rds"), emit: signatures_nmfOut_rds
        tuple val(meta), path("*_plot_all.rds"), emit: signatures_plot_rds
        tuple val(meta), path("*_plot_all.pdf"), emit: signatures_plot_pdf
        path "versions.yml", emit: versions

    script:
    template "main_script.R"
}
