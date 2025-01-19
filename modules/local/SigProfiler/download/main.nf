process DOWNLOAD_GENOME_SIGPROFILER {
    container = 'docker://katiad/sigprofiler:version1.0'

    input:
        val(reference_genome) // reference_genome : genome -> for example: GRCh37

    output:
        path("*"), emit: genome_sigprofiler
        path "versions.yml", emit: versions


    script:
    """
    SigProfilerMatrixGenerator install $reference_genome -v .


    VERSION=\$(pip show SigProfilerMatrixGenerator | grep Version | awk '{print \$NF}')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SigProfilerMatrixGenerator: \$VERSION
    END_VERSIONS

    """


}
