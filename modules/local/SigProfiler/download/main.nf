process DOWNLOAD_GENOME_SIGPROFILER {
    container = 'docker://katiad/sigprofiler:version1.0'

    input:
        val(reference_genome) // reference_genome : genome -> for example: GRCh37

    output:
        path("*"), emit: genome_sigprofiler


    script:
    """
    SigProfilerMatrixGenerator install $reference_genome -v .

    """

    stub:
    """
    echo "${task.process}:" > versions.yml
    echo 'SigProfilerMatrixGenerator:v1.2.29' >> versions.yml
    """

}
