workflow ANNOTATION_CACHE_INITIALISATION {
    take:
    vep_cache
    vep_species
    vep_cache_version
    vep_genome
    vep_custom_args

    main:
    def vep_annotation_cache_key = (vep_cache == "s3://annotation-cache/vep_cache/") ? "${vep_cache_version}_${vep_genome}/" : ""
    def vep_species_suffix = vep_custom_args.contains("--merged") ? '_merged' : (vep_custom_args.contains("--refseq") ? '_refseq' : '')
    def vep_cache_dir = "${vep_annotation_cache_key}${vep_species}${vep_species_suffix}/${vep_cache_version}_${vep_genome}"
    def vep_cache_path_full = file("$vep_cache/$vep_cache_dir", type: 'dir')
    
    if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
        if (vep_cache == "s3://annotation-cache/vep_cache/") {
            error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
        } else {
            error("Path provided with VEP cache is invalid.\nMake sure there is a directory named ${vep_cache_dir} in ${vep_cache}./n")
        }
    }
    ensemblvep_cache = Channel.fromPath(file("${vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()

    emit:
    ensemblvep_cache
}
