include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'

workflow DOWNLOAD_CACHE_VEP {
    take:
        ensemblvep_info

    main:
        ENSEMBLVEP_DOWNLOAD(ensemblvep_info)
    
    emit:
        ensemblvep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect()  // channel: [ meta, cache ]
}