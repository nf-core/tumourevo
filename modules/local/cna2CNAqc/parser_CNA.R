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
    purity = solutions$cellularity[2]
    ploidy = solutions$ploidy.estimate[2]
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
    purity = solutions$AberrantCellFraction
    ploidy = solutions$Ploidy
    return(list(segments = segments, purity = purity, ploidy = ploidy))
}
