#!/usr/bin/env Rscript

parse_args = function(x){
    args_list = unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals = lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))
    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals = lapply(args_vals, function(z){ length(z) = 2; z})
    parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    purple_max_cn_deviation = 0.25
)
args_opt = parse_args('$task.ext.args')
for (ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]

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
    purity = solutions[["cellularity"]][1]
    ploidy = solutions[["ploidy.estimate"]][1]
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

parse_facets = function(segments_file, extra_file){
    data = vcfR::read.vcfR(segments_file, verbose = FALSE)
    meta = data@meta
    purity_line = meta[grep("^##purity=", meta)]
    ploidy_line = meta[grep("^##ploidy=", meta)]

    purity_d = if (length(purity_line)) as.numeric(sub("^##purity=", "", purity_line[1]))
    ploidy_d = if (length(ploidy_line)) as.numeric(sub("^##ploidy=", "", ploidy_line[1]))

    fix = as.data.frame(vcfR::getFIX(data), stringsAsFactors = FALSE)
    # INFO fields
    endpos = as.integer(vcfR::extract.info(data, "END"))
    tcn_em = as.integer(vcfR::extract.info(data, "TCN_EM"))
    lcn_em = as.integer(vcfR::extract.info(data, "LCN_EM"))

    seg = data.frame(
        chr      = fix[["CHROM"]],
        from = as.integer(fix[["POS"]]),
        to   = endpos,
        Major   = ifelse(is.na(tcn_em) | is.na(lcn_em), NA_integer_, tcn_em - lcn_em),
        minor   = lcn_em
    )

    seg = seg[
        !is.na(seg[["to"]]) &
        !is.na(seg[["Major"]]) &
        !is.na(seg[["minor"]]),
    ]

    return(list(segments = seg, purity = purity_d, ploidy = ploidy_d))
}

parse_PURPLE = function(segments_file, extra_file, max_deviation = 0.25){
    # PURPLE reports continuous allele copy numbers; CNAqc needs integer clonal karyotypes.
    # Segments where an allele copy number is far from an integer (subclonal or noisy copy
    # number) cannot be represented as a clonal karyotype: rounding them assigns the wrong
    # copy number to their mutations downstream, so they are removed instead of rounded.
    # The default threshold of 0.25 uses the same numerical cutoff reported by PURPLE's
    # polyclonalProportion. This downstream rule is deliberately conservative: it is
    # applied to either allele and does not reproduce PURPLE's confidence weighting.
    if (length(max_deviation) != 1 || is.na(max_deviation) ||
        max_deviation < 0 || max_deviation > 0.5) {
        stop("PURPLE maximum copy-number deviation must be a number between 0 and 0.5.")
    }

    segments = readr::read_tsv(segments_file, col_types = readr::cols()) %>%
        dplyr::mutate(
            chr = chromosome,
            from = start,
            to = end,
            Major = as.integer(round(majorAlleleCopyNumber)),
            minor = as.integer(round(minorAlleleCopyNumber)),
            deviation = pmax(abs(majorAlleleCopyNumber - Major), abs(minorAlleleCopyNumber - minor))) %>%
        dplyr::filter(!is.na(Major), !is.na(minor), Major >= 0, minor >= 0)

    if (nrow(segments) == 0) {
        stop("No valid non-negative PURPLE copy-number segments were found.")
    }

    off_integer = segments[["deviation"]] > max_deviation
    segment_bp = segments[["to"]] - segments[["from"]] + 1
    message(sprintf(
        "PURPLE: removing %d of %d segments (%.1f%% of the segmented genome) with an allele copy number more than %.2f from an integer",
        sum(off_integer), length(off_integer), 100 * sum(segment_bp[off_integer]) / sum(segment_bp), max_deviation))
    if (all(off_integer)) {
        stop("All PURPLE segments have an allele copy number more than ", max_deviation,
             " from an integer; no clonal copy number left for CNAqc (see --purple_max_cn_deviation).")
    }
    removed_segments = segments[off_integer, ] %>%
        dplyr::select(chr, from, to, majorAlleleCopyNumber, minorAlleleCopyNumber, deviation)
    segments = segments[!off_integer, ] %>%
        dplyr::select(chr, from, to, Major, minor)

    solutions = readr::read_tsv(extra_file, col_types = readr::cols())
    purity = solutions[["purity"]][1]
    ploidy = solutions[["ploidy"]][1]

    # The status column of purple.purity.tsv holds PURPLE's fit method: NORMAL,
    # HIGHLY_DIPLOID, SOMATIC or NO_TUMOR. PURPLE's QC verdict (QCStatus: PASS, WARN_*,
    # FAIL_*) is in <tumour>.purple.qc, which is not an input; check it separately.
    fit_method = solutions[["status"]][1]
    if (!is.null(fit_method) && !is.na(fit_method)) {
        if (fit_method == "NO_TUMOR") {
            stop("PURPLE fit method is NO_TUMOR: no tumour was detected, so its copy number cannot be used.")
        }
        if (fit_method != "NORMAL") {
            warning("PURPLE fit method is ", fit_method, "; purity, ploidy and copy number may be unreliable.")
        }
    }
    if (!is.na(purity) && purity < 0.2) {
        warning("PURPLE purity is ", purity, ", below 0.2; copy number may be unreliable.")
    }

    return(list(segments = segments, purity = purity, ploidy = ploidy, removed_segments = removed_segments))
}

if ("$meta.cna_caller" == 'sequenza'){
CNA = parse_Sequenza(segments = "$cna_segs", extra = "$cna_extra")

} else if ("$meta.cna_caller" == 'ASCAT'){
CNA = parse_ASCAT(segments = "$cna_segs", extra = "$cna_extra")

} else if ("$meta.cna_caller" == 'Battenberg'){
CNA = parse_Battenberg(segments = "$cna_segs", extra = "$cna_extra")

} else if ("$meta.cna_caller" == 'facets'){
CNA = parse_facets(segments = "$cna_segs", extra = "$cna_extra")

} else if ("$meta.cna_caller" == 'PURPLE'){
CNA = parse_PURPLE(segments = "$cna_segs", extra = "$cna_extra",
                   max_deviation = as.numeric(opt[["purple_max_cn_deviation"]]))

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
