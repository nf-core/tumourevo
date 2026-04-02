process CNAQC2TSV {
    tag "$meta.id"
    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/913112a2d9295d35fe97874caf5f500df5a98d6fef1cb5861fd64caa0223a047/data':
        'community.wave.seqera.io/library/r-cnaqc_r-cli_r-dplyr_r-readr_pruned:0fc82bfd06afe6dc' }"

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
