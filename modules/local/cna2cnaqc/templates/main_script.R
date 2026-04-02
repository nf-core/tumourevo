#!/usr/bin/env Rscript

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)

library(dplyr)
library(readr)

parse_Sequenza = function(segments_file, extra_file){
    # Extract the segments information
    segments = readr::read_tsv(segments_file, col_types = readr::cols()) %>%
        dplyr::rename(
            chr = chromosome,
            from = start.pos,
            to = end.pos,
            Major = A,
            minor = B) %>%
        dplyr::select(chr, from, to, Major, minor, dplyr::everything())

    solutions = readr::read_tsv(extra_file, col_types = readr::cols())
    purity = solutions[["cellularity"]][2]
    ploidy = solutions[["ploidy.estimate"]][2]
    return(list(segments = segments, purity = purity, ploidy = ploidy))
}

parse_ASCAT = function(segments_file, extra_file){
    # Extract the segments information
    segments = readr::read_tsv(segments_file, col_types = readr::cols()) %>%
        dplyr::mutate(chr = paste0("chr",chr)) %>%
        dplyr::rename(
            from = startpos,
            to = endpos,
            Major = nMajor,
            minor = nMinor) %>%
        dplyr::select(chr, from, to, Major, minor)

    solutions = readr::read_tsv(extra_file, col_types = readr::cols())
    purity = solutions[["AberrantCellFraction"]]
    ploidy = solutions[["Ploidy"]]
    return(list(segments = segments, purity = purity, ploidy = ploidy))
}

parse_Battenberg = function(segments_file, extra_file){
    # Extract the extra information
    purity_ploidy = read.table(extra_file, sep = '\t', header = T) %>%
        dplyr::filter(is.best == TRUE) %>%
        dplyr::rename(purity = rho) %>%
        dplyr::select(purity, ploidy)

    # Extract the segments information
    segments = readr::read_tsv(segments_file, col_types = readr::cols()) %>%
        dplyr::select(chr, startpos, endpos, nMaj1_A, nMin1_A, nMaj2_A, nMin2_A,frac1_A) %>%
        dplyr::rename(from = startpos, to = endpos) %>%
        dplyr::rename(Major = nMaj1_A, minor = nMin1_A, Major_2 = nMaj2_A, minor_2 = nMin2_A, CCF = frac1_A)

    return(list(segments = segments, purity = purity_ploidy[['purity']], ploidy = purity_ploidy[['ploidy']]))
}

if ("$meta.cna_caller" == 'sequenza'){
CNA = parse_Sequenza(segments = "$cna_segs", extra = "$cna_extra")

} else if ("$meta.cna_caller" == 'ASCAT'){
CNA = parse_ASCAT(segments = "$cna_segs", extra = "$cna_extra")

} else if ("$meta.cna_caller" == 'Battenberg'){
CNA = parse_Battenberg(segments = "$cna_segs", extra = "$cna_extra")

} else {
stop('Copy Number Caller not supported.')
}

saveRDS(object = CNA, file = paste0(opt[["prefix"]], "_cna.rds"))

# version export
f = file("versions.yml","w")
readr_version = sessionInfo()\$otherPkgs\$readr\$Version
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    readr:", readr_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
close(f)
