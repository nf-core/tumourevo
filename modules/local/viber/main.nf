process VIBER {
    tag "$meta.id"
    label "process_single"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://elenabuscaroli/viber:version0.1' :
        'docker.io/elenabuscaroli/viber:version0.1' }"

    input:
        tuple val(meta), path(rds_join), val(tumour_samples) //rds from either JOIN_CNAQC or JOIN_FIT, should be always grouped

    output:
        tuple val(meta), path("*_viber_best_st_fit.rds"), emit: viber_rds
        tuple val(meta), path("*_viber_best_st_heuristic_fit.rds"), emit: viber_heuristic_rds
        tuple val(meta), path("*_${plot1}"), emit: viber_plots_rds
        tuple val(meta), path("*_${plot2}"), emit: viber_heuristic_plots_rds
        tuple val(meta), path("*_REPORT_plots_viber.rds"), emit: viber_report_rds
        tuple val(meta), path("*_REPORT_plots_viber.pdf"), emit: viber_report_pdf
        tuple val(meta), path("*_REPORT_plots_viber.png"), emit: viber_report_png
        path "versions.yml", emit: versions

    script:
    def n_samples = tumour_samples.size()
    if (n_samples==1) {
        plot1 = "viber_best_st_mixing_plots.rds"
        plot2 = "viber_best_st_heuristic_mixing_plots.rds"
    } else {
        plot1 = "viber_best_st_fit_plots.rds"
        plot2 = "viber_best_st_heuristic_fit_plots.rds"
    }

    template "main_script.R"
}
