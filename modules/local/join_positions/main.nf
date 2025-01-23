process JOIN_POSITIONS {
    tag "$meta.id"
    label "process_single"
    container = 'docker://lvaleriani/cnaqc:version1.0'

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
    tb = vcfR::vcfR2tidy(pileup)

    gt_field = tb\$gt %>%
            tidyr::separate(gt_AD, sep = ',', into = c("NR", "NV")) %>%
            dplyr::rename(DP = gt_DP) %>%
            dplyr::mutate(
                NR = as.numeric(NR),
                NV = as.numeric(NV),
                VAF = NV/DP) %>%
            dplyr::mutate(VAF = ifelse(is.na(VAF),0,VAF)) %>%
            dplyr::rename(sample = Indiv)

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
                dplyr::select(chr, from, to, ref, NV, DP, VAF, dplyr::everything(), -alt, -ChromKey) %>%
                dplyr::mutate(from = from+1, to = to+1)  %>%
                dplyr::mutate(NV = 0)

    bind = inner_join(pileup_mutations, all_positions, by = join_by(chr, from, to, ref))

    vcf[["$meta.tumour_sample"]]\$mutations = dplyr::bind_rows(mutations, bind)
    saveRDS(file = paste0("$prefix", "_pileup_VCF.rds"), object = vcf)


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
}
