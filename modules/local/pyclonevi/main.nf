process PYCLONEVI {
    tag "$meta.id"
    label "process_low"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyclone-vi%3A0.1.3--pyhca03a8a_0' :
        'docker.io/blcdsdockerregistry/pyclone-vi:0.1.2' }"

    input:
        tuple val(meta), path(rds_join), val(tumour_samples)

    output:
        tuple val(meta), path("*_cluster_table.csv"), emit: ctree_input
        tuple val(meta), path("*.tsv"), emit: pyclone_input
        tuple val(meta), path("*_all_fits.h5"), emit: pyclone_all_fits
        tuple val(meta), path("*_best_fit.txt"), emit: pyclone_best_fit
        path "versions.yml", emit: versions

    script:
    template "main_script.py"
}
