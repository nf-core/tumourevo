process JOIN_POSITIONS {
    tag "$meta.id"
    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/913112a2d9295d35fe97874caf5f500df5a98d6fef1cb5861fd64caa0223a047/data':
        'community.wave.seqera.io/library/r-cnaqc_r-cli_r-dplyr_r-readr_pruned:0fc82bfd06afe6dc' }"

    input:
    tuple val(meta), path(rds), path(vcf_pileup), path(positions)

    output:
    tuple val(meta), path("*.rds"),     emit: rds
    path "versions.yml",                emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(tidyr)
    library(dplyr)
    library(vcfR)

    all_positions = readRDS("$positions")

    vcf = readRDS("$rds")
    mutations = vcf[["$meta.tumour_sample"]]\$mutations

    pileup = vcfR::read.vcfR("$vcf_pileup")
    if (nrow(pileup@fix) == 0) {

    message("VCF is empty: no variants found, original vcf is returned.")
    saveRDS(file = paste0("$prefix", "_pileup_VCF.rds"), object = vcf)
    } else {

    tb = vcfR::vcfR2tidy(pileup)
    gt_field = tb\$gt %>%
            tidyr::separate(gt_AD, sep = ',', into = c("NR", "NV")) %>%
            dplyr::rename(DP = gt_DP) %>%
            dplyr::mutate(
                NR = as.numeric(NR),
                NV = as.numeric(NV),
                VAF = NV/DP) %>%
            dplyr::mutate(VAF = ifelse(is.na(VAF),0,VAF)) %>%
	    dplyr::mutate(sample = "$meta.tumour_sample") %>%
            dplyr::select(!Indiv)

    fix_field = tb\$fix %>%
            dplyr::rename(
                chr = CHROM,
                from = POS,
                ref = REF,
                alt = ALT
            ) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                from = as.numeric(from) - 1,
                to = from + 1
            ) %>%
            dplyr::ungroup() %>%
            dplyr::select(chr, from, to, ref, alt, AD)

    if(nrow(fix_field) != nrow(gt_field))
        stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")

    pileup_mutations = dplyr::bind_cols(fix_field, gt_field)  %>%
                dplyr::select(chr, from, to, ref, NV, DP, dplyr::everything(), -alt, -ChromKey) %>%
                dplyr::mutate(from = from+1, to = to+1)  %>%
                dplyr::mutate(NV = 0) %>%
                dplyr::mutate(VAF = NV/DP)

    bind = inner_join(pileup_mutations, all_positions, by = join_by(chr, from, to, ref))

    vcf[["$meta.tumour_sample"]]\$mutations = dplyr::bind_rows(mutations, bind)
    saveRDS(file = paste0("$prefix", "_pileup_VCF.rds"), object = vcf)
}

    # version export
    f <- file("versions.yml","w")
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    tidyr_version <- sessionInfo()\$otherPkgs\$tidyr\$Version
    vcfR_version <- sessionInfo()\$otherPkgs\$vcfR\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    writeLines(paste("    tidyr:", tidyr_version), f)
    writeLines(paste("    vcfR:", vcfR_version), f)
    close(f)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pileup_VCF.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        tidyr: \$(Rscript -e "cat(as.character(packageVersion('tidyr')))")
        vcfR: \$(Rscript -e "cat(as.character(packageVersion('vcfR')))")
    END_VERSIONS
    """
}
