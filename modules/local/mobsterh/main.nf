process MOBSTERh {
    tag "$meta.id"
    label "process_single"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://elenabuscaroli/mobster:version1.0' :
        'docker.io/elenabuscaroli/mobster:version1.0' }"

    input:
        tuple val(meta), path(rds_join) // rds from JOIN_CNAQC

    output:
        tuple val(meta), path("*_mobsterh_st_fit.rds"), emit: mobster_rds
        tuple val(meta), path("*_mobsterh_st_best_fit.rds"), emit: mobster_best_rds
        tuple val(meta), path("*_plots.rds"), emit: mobster_plots_rds
        tuple val(meta), path("*_REPORT_plots_mobster.rds"), emit: mobster_report_rds
        tuple val(meta), path("*_REPORT_plots_mobster.pdf"), emit: mobster_report_pdf
        tuple val(meta), path("*_REPORT_plots_mobster.png"), emit: mobster_report_png
        path "versions.yml", emit: versions

    script:
    template "main_script.R"
}
