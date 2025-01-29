#!/usr/bin/env Rscript

parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
    # giving errors when we have lists like c(xxx, xxx) since it will separate it
    # args_list = unlist(strsplit(x, ', ')[[1]])
    args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
    # args_vals = lapply(args_list, function(x) strsplit(x, split=":")[[1]])
    args_vals = lapply(args_list, function(x) {
        x_splt = strsplit(x, split=":")[[1]]
        c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
    })

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

    parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)
args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]

print("\n\n\n")
print(opt)
print("\n\n\n")
print("$task.ext.args")
print("\n\n\n")

# Auxiliriaty functions ####

parse_FreeBayes = function(vcf, tumour_id, normal_id, filter_mutations = FALSE) {
    tb = vcfR::vcfR2tidy(vcf)

    gt_field = tb[["gt"]] %>%
        tidyr::separate(gt_AD, sep = ",", into = c("NR", "NV")) %>%
        dplyr::mutate(
            NR = as.numeric(NR),
            NV = as.numeric(NV),
            DP = NV + NR,
            VAF = NV/DP) %>%
        dplyr::rename(sample = Indiv)

    samples_list = gt_field[["sample"]] %>% unique

    if ("CSQ" %in% tb[["meta"]][["ID"]]){
    # VEP specific field extraction
    # Take CSQ field names and split by |

        vep_field = tb[["meta"]] %>%
                dplyr::filter(ID == "CSQ") %>%
                dplyr::select(Description) %>%
                dplyr::pull()

        vep_field = strsplit(vep_field, split = "|", fixed = TRUE)[[1]]
        print(vep_field)
        vep_field = vep_field[2:length(vep_field)-1]

    # Tranform the fix field by splittig the CSQ and select the columns needed
        fix_field = tb[["fix"]] %>%
            dplyr::rename(
                chr = CHROM,
                from = POS,
                ref = REF,
                alt = ALT) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                from = as.numeric(from),
                to = from + nchar(alt)) %>%
            dplyr::ungroup() %>%
            dplyr::select(chr, from, to, ref, alt, CSQ, dplyr::everything(),  -ChromKey, -DP) %>%
            tidyr::separate(CSQ, vep_field, sep = "|") %>%
            dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything())
        print(fix_field)

    } else {
        fix_field = tb[["fix"]] %>%
            dplyr::rename(
            chr = CHROM,
            from = POS,
            ref = REF,
            alt = ALT
            ) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
            from = as.numeric(from),
            to = from + nchar(alt)) %>%
            dplyr::ungroup() %>%
            dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey, -DP)
    }

    # if have to filter mutations
    if (filter_mutations){
        filter = c('PASS')
    } else {
        filter = fix_field[["FILTER"]] %>% unique()
    }

    calls = lapply(
        samples_list,
        function(s){
        gt_field_s = gt_field %>% dplyr::filter(sample == s)

        if(nrow(fix_field) != nrow(gt_field_s))
            stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")

        fits = list()
        fits[["sample"]] = s
        fits[["mutations"]] = dplyr::bind_cols(fix_field, gt_field_s) %>%
            dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything()) %>%
            dplyr::filter(FILTER %in% filter)
        fits
        })

    names(calls) = samples_list
    samples = c(tumour_id, normal_id)
    calls = calls[samples]
    return(calls)
}


parse_Mutect = function(vcf, tumour_id, normal_id, filter_mutations = FALSE){
    # Transform vcf to tidy
    tb = vcfR::vcfR2tidy(vcf)

    # Extract gt field and obtain coverage (DP) and variant allele frequency (VAF) fields
    gt_field = tb[["gt"]] %>%
        tidyr::separate(gt_AD, sep = ",", into = c("NR", "NV")) %>%
        dplyr::mutate(
            NR = as.numeric(NR),
            NV = as.numeric(NV),
            DP = NV + NR,
            VAF = NV/DP) %>%
        dplyr::rename(sample = Indiv)

    # Extract sample names
    samples_list = gt_field[["sample"]] %>% unique

    # check if VCF is annotated with VEP
    if ("CSQ" %in% tb[["meta"]][["ID"]]){
    # VEP specific field extraction
    # Take CSQ field names and split by |

        vep_field = tb[["meta"]] %>%
                    dplyr::filter(ID == "CSQ") %>%
                    dplyr::select(Description) %>%
                    dplyr::pull()

        vep_field = strsplit(vep_field, split = "|", fixed = TRUE)[[1]]
        print(vep_field)
        vep_field = vep_field[2:length(vep_field)-1]

    # Tranform the fix field by splittig the CSQ and select the columns needed
        fix_field = tb[["fix"]] %>%
            dplyr::rename(
                chr = CHROM,
                from = POS,
                ref = REF,
                alt = ALT) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                from = as.numeric(from),
                to = from + nchar(alt)) %>%
            dplyr::ungroup() %>%
            dplyr::select(chr, from, to, ref, alt, CSQ, dplyr::everything()) %>%
            tidyr::separate(CSQ, vep_field, sep = "|") %>%
            dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything(), -DP) #can add other thing, CSQ, HGSP
        print(fix_field)

    } else {
        # Take from fix field some columns
        fix_field = tb[["fix"]] %>%
            dplyr::rename(
                chr = CHROM,
                from = POS,
                ref = REF,
                alt = ALT) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                from = as.numeric(from),
                to = from + nchar(alt)) %>%
            dplyr::ungroup() %>%
            dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey, -DP) #-DP
    }

    if (filter_mutations){
        filter = c('PASS')
    } else {
        filter = fix_field[["FILTER"]] %>% unique()
    }

    # For each sample create the table of mutations
    calls = lapply(
            samples_list,
            function(s){
            gt_field_s = gt_field %>% dplyr::filter(sample == s)

            if(nrow(fix_field) != nrow(gt_field_s))
                stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")

            fits = list()
            fits[["sample"]] = s
            fits[["mutations"]] = dplyr::bind_cols(fix_field, gt_field_s) %>%
                dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything()) %>%
                dplyr::filter(FILTER %in% filter)
            fits
            })

    names(calls) = samples_list
    samples = c(tumour_id, normal_id)
    calls = calls[samples]
    return(calls)
}

retrieve_ref_alt = function(row){
    ref = row[["ref"]]
    alt = row[["alt"]]
    if (ref == 'A'){NR=as.integer(strsplit(row[["gt_AU"]], split=',')[[1]][1])}
    if (ref == 'T'){NR=as.integer(strsplit(row[["gt_TU"]], split=',')[[1]][1])}
    if (ref == 'G'){NR=as.integer(strsplit(row[["gt_GU"]], split=',')[[1]][1])}
    if (ref == 'C'){NR=as.integer(strsplit(row[["gt_CU"]], split=',')[[1]][1])}

    if (alt == 'A'){NV=as.integer(strsplit(row[["gt_AU"]], split=',')[[1]][1])}
    if (alt == 'T'){NV=as.integer(strsplit(row[["gt_TU"]], split=',')[[1]][1])}
    if (alt == 'G'){NV=as.integer(strsplit(row[["gt_GU"]], split=',')[[1]][1])}
    if (alt == 'C'){NV=as.integer(strsplit(row[["gt_GU"]], split=',')[[1]][1])}

    ref_alt = paste0(NR, ',', NV)
    ref_alt
}

parse_Strelka = function(vcf, tumour_id, normal_id, filter_mutations = FALSE){
    tb = vcfR::vcfR2tidy(vcf)
    gt_field = tb[["gt"]] %>% rename(sample = Indiv)
    samples_list = gt_field[["sample"]] %>% unique

    if ("CSQ" %in% tb[["meta"]][["ID"]]){
    # VEP specific field extraction
    # Take CSQ field names and split by |

        vep_field = tb[["meta"]] %>%
                dplyr::filter(ID == "CSQ") %>%
                dplyr::select(Description) %>%
                dplyr::pull()

    vep_field = strsplit(vep_field, split = "|", fixed = TRUE)[[1]]
    print(vep_field)
    vep_field = vep_field[2:length(vep_field)-1]

    # Tranform the fix field by splittig the CSQ and select the columns needed
    fix_field = tb[["fix"]] %>%
        dplyr::rename(
            chr = CHROM,
            from = POS,
            ref = REF,
            alt = ALT) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            from = as.numeric(from),
            to = from + nchar(alt)) %>%
        dplyr::ungroup() %>%
        dplyr::select(chr, from, to, ref, alt, CSQ, dplyr::everything(), -ChromKey) %>%
        tidyr::separate(CSQ, vep_field, sep = "|") %>%
        dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything()) #can add other thing, CSQ, HGSP
    print(fix_field)

    } else {
        # Take from fix field some columns
        fix_field = tb[["fix"]] %>%
            dplyr::rename(
                chr = CHROM,
                from = POS,
                ref = REF,
                alt = ALT) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                from = as.numeric(from),
                to = from + nchar(alt)) %>%
            dplyr::ungroup() %>%
            dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey) #-DP
    }

    if (filter_mutations){
        filter = c('PASS')
    } else {
        filter = fix_field[["FILTER"]] %>% unique()
    }

    calls = lapply(
        samples_list,
        function(s){
            gt_field_s = gt_field %>% dplyr::filter(sample == s)

            fits = list()
            fits[["sample"]] = s
            mutations = dplyr::bind_cols(fix_field, gt_field_s)
            ref_alt = lapply(1:nrow(mutations), function(r){
                retrieve_ref_alt(mutations[r,])
        })

        ref_alt = ref_alt %>% unlist()
        mutations[["ref_alt"]] = ref_alt
        mutations = mutations %>%
            tidyr::separate(ref_alt, into = c('NR', 'NV')) %>%
            dplyr::mutate(NR = as.integer(NR),
                        NV = as.integer(NV)) %>%
            dplyr::mutate(DP = NR+NV) %>%
            dplyr::mutate(VAF = NV/DP) %>%
            dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, everything()) %>%
            dplyr::filter(FILTER %in% filter)

        fits[["mutations"]] = mutations
        fits
        })

    names(calls) = samples_list
    samples = c(tumour_id, normal_id)
    calls = calls[samples]
    return(calls)
}


parse_Platypus = function(vcf, tumour_id, normal_id, filter_mutations = FALSE){
    tb = vcfR::vcfR2tidy(vcf)

    gt_field = tb[["gt"]] %>%
        dplyr::mutate(
            DP = as.numeric(gt_NR),
            NV = as.numeric(gt_NV),
            VAF = NV/DP) %>%
        dplyr::rename(sample = Indiv)

    # Extract each sample
    samples_list = gt_field[["sample"]] %>% unique

    # check if VCF is annotated with VEP
    if ("CSQ" %in% tb[["meta"]][["ID"]]){
    # VEP specific field extraction
    # Take CSQ field names and split by |

        vep_field = tb[["meta"]] %>%
                dplyr::filter(ID == "CSQ") %>%
                dplyr::select(Description) %>%
                dplyr::pull()

    vep_field = strsplit(vep_field, split = "|", fixed = TRUE)[[1]]
    print(vep_field)
    vep_field = vep_field[2:length(vep_field)-1]

    # Tranform the fix field by splittig the CSQ and select the columns needed
    fix_field = tb[["fix"]] %>%
        dplyr::rename(
            chr = CHROM,
            from = POS,
            ref = REF,
            alt = ALT) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            from = as.numeric(from),
            to = from + nchar(alt)) %>%
        dplyr::ungroup() %>%
        dplyr::select(chr, from, to, ref, alt, CSQ, dplyr::everything(),  -ChromKey) %>%
        tidyr::separate(CSQ, vep_field, sep = "|") %>%
        dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything())
    print(fix_field)

    } else {
        fix_field = tb[["fix"]] %>%
            dplyr::rename(
            chr = CHROM,
            from = POS,
            ref = REF,
            alt = ALT) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
            from = as.numeric(from),
            to = from + nchar(alt)) %>%
            dplyr::ungroup() %>%
            dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey)
    }

    if (filter_mutations){
        filter = c('PASS')
    } else {
        filter = fix_field[["FILTER"]] %>% unique()
    }

    calls = lapply(
        samples_list,
        function(s){
        gt_field_s = gt_field %>%
            dplyr::filter(sample == s)

        if(nrow(fix_field) != nrow(gt_field_s))
            stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")

        fits = list()
        fits[["sample"]] = s
        fits[["mutations"]] = dplyr::bind_cols(fix_field, gt_field_s) %>%
            dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything()) %>%
            dplyr::filter(FILTER %in% filter)
        fits
        })

    names(calls) = samples_list
    samples = c(tumour_id, normal_id)
    calls = calls[samples]
    return(calls)
}


library(dplyr)
library(tidyr)
library(vcfR)


# Read vcf file
vcf = vcfR::read.vcfR("$vcf")

# Check from which caller the .vcf has been produced
source = vcfR::queryMETA(vcf, element = 'source')[[1]]

if (TRUE %in% grepl(pattern = 'Mutect', x = source)){
    calls = parse_Mutect(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample", filter_mutations = as.logical(opt[["filter_mutations"]]))

} else if (TRUE %in% grepl(pattern = 'strelka', x = source)){
    calls = parse_Strelka(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample", filter_mutations = as.logical(opt[["filter_mutations"]]))

} else if (TRUE %in% grepl(pattern = 'Platypus', x = source)){
    calls = parse_Platypus(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample", filter_mutations = as.logical(opt[["filter_mutations"]]))

} else if (TRUE %in% grepl(pattern = 'freeBayes', x = source)){
    calls = parse_FreeBayes(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample", filter_mutations = as.logical(opt[["filter_mutations"]]))

} else {
    stop('Variant Caller not supported.')
}

saveRDS(object = calls, file = paste0(opt[["prefix"]], "_snv.rds"))

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
