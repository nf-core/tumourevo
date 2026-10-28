#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/tumourevo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/tumourevo
    Website: https://nf-co.re/tumourevo
    Slack  : https://nfcore.slack.com/channels/tumourevo
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TUMOUREVO } from './workflows/tumourevo'
include { samplesheetToList } from 'plugin/nf-schema'

include { ANNOTATION_CACHE_INITIALISATION  } from './subworkflows/local/annotation_cache_initialisation'
include { ENSEMBLVEP_DOWNLOAD } from './modules/nf-core/ensemblvep/download/main'
include { UNTAR } from './modules/nf-core/untar'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_TUMOUREVO {

    take:
    input

    main:
    ch_versions = Channel.empty()
    fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
    drivers_table =  params.drivers_table ? Channel.fromPath(params.drivers_table).map{ it -> [ it ] }.collect() : Channel.empty()

    if (params.download_cache_vep) {
        ensemblvep_info = Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
        ENSEMBLVEP_DOWNLOAD(ensemblvep_info)
        vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect().map{ meta, cache -> [ cache ] }

    } else {
        if (params.vep_cache.endsWith("tar.gz")){
            def vep_tar = params.vep_cache ? Channel.fromPath(params.vep_cache).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
            UNTAR(vep_tar) //[meta, path]
            vep_cache = UNTAR.out.untar.map{meta, cache_path ->
                [cache_path]}

        } else {

            ANNOTATION_CACHE_INITIALISATION(
                params.vep_cache,
                params.vep_species,
                params.vep_cache_version,
                params.vep_genome,
                params.vep_custom_args)

            vep_cache  = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache
        }
    }

    TUMOUREVO (
        input,
        fasta,
        drivers_table,
        vep_cache,
        ch_versions
    )
    ch_versions = ch_versions.mix(TUMOUREVO.out.versions)

    emit:
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    def input = params.input
      ? Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))
      : Channel.empty()

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_TUMOUREVO(
        input
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
