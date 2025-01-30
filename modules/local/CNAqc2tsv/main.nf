process RDS_PROCESSING {
    tag "$meta.id"
    label "process_single"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/cnaqc:version1.0' :
        'docker.io/lvaleriani/cnaqc:version1.0' }"

    input:
    tuple val(meta), path(rds_join), val(tumour_samples)

    output:
    tuple val(meta), path("*_joint_table.tsv"), val(tumour_samples),    emit: tsv
    path "versions.yml",                                                emit: versions

    script:
    template "main_script.R"
}
