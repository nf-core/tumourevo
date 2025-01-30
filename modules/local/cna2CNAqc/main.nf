process CNA_PROCESSING {
    tag "$meta.id"
    label "process_single"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/cnaqc:version1.0' :
        'docker.io/lvaleriani/cnaqc:version1.0' }"

    input:
    tuple val(meta), path(cna_segs), path(cna_extra)

    output:
    tuple val(meta), path("*.rds"),                            emit: rds
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"
}
