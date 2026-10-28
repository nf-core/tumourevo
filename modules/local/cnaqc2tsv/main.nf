process CNAQC2TSV {
    tag "$meta.id"
    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/63f20f0f8ea78a7aecb1baa9aa35c03cdff9e18d19181d34561ed6d7e6376a0c/data' :
        'community.wave.seqera.io/library/r-cnaqc_r-cli_r-dplyr_r-readr_pruned:aebdc556849a14c7' }"

    input:
    tuple val(meta), path(rds_join), val(tumour_samples)

    output:
    tuple val(meta), path("*_joint_table.tsv"), val(tumour_samples), emit: tsv
    path "versions.yml",                                             emit: versions

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_joint_table.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnaqc: \$(Rscript -e "cat(as.character(packageVersion('CNAqc')))")
        readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
