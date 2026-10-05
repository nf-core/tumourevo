#!/usr/bin/env Rscript

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)

parse_Mutect = function(vcf, tumour_id, normal_id){
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
        dplyr::select(chr, from, to, ref, alt, dplyr::everything(), -ChromKey, -DP)

    # Extract sample names
    samples_list = gt_field[["sample"]] %>% unique

    calls = lapply(
        samples_list,
        function(s){
            gt_field_s = gt_field %>% dplyr::filter(sample == s)

            if(nrow(fix_field) != nrow(gt_field_s))
                stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")

            fits = list()
            fits[["sample"]] = s
            fits[["mutations"]] = dplyr::bind_cols(fix_field, gt_field_s) %>%
                dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything())
        fits
    })

    names(calls) = samples_list
    if (normal_id == '[]'){
      samples = c(tumour_id)
    } else{
      samples = c(tumour_id, normal_id)
    }
    calls = calls[samples]

    # check if VCF is annotated with VEP
    if ("CSQ" %in% tb[["meta"]][["ID"]]){
    # VEP specific field extraction
    # Take CSQ field names and split by |

        vep_field = tb[['meta']] %>%
            dplyr::filter(ID == "CSQ") %>%
            dplyr::select(Description) %>%
            dplyr::pull()

        tmp_vep_field = strsplit(vep_field, split = "|", fixed = TRUE) %>% unlist()
        vep_field = tmp_vep_field[1:length(tmp_vep_field)-1]

        # Tranform the fix field by splittig the CSQ and select the columns needed
        calls[[tumour_id]][['mutations']] = calls[[tumour_id]][['mutations']] %>%
            dplyr::mutate(CSQ = strsplit(CSQ, ",")) %>%
            tidyr::unnest(CSQ) %>%
            tidyr::separate(CSQ, vep_field, sep = "\\\\|") %>%
            dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything())  #can add other thing, CSQ, HGSP

        if (normal_id != '[]'){
          calls[[normal_id]][['mutations']] = calls[[normal_id]][['mutations']] %>% dplyr::select(-CSQ) %>% dplyr::distinct()
        }
    }
    return(calls)
}

retrieve_ref_alt = function(row){
    ref = row[["ref"]]
    alt = row[["alt"]]
    if (ref == 'A'){NR=as.integer(strsplit(row[["gt_AU"]], split=',')[[1]][1])}
    if (ref == 'T'){NR=as.integer(strsplit(row[["gt_TU"]], split=',')[[1]][1])}
    if (ref == 'G'){NR=as.integer(strsplit(row[["gt_GU"]], split=',')[[1]][1])}
    if (ref == 'C'){NR=as.integer(strsplit(row[["gt_CU"]], split=',')[[1]][1])}
    if (nchar(ref) > 1 | nchar(alt) > 1 ){NR=as.integer(strsplit(row[["gt_TAR"]], split=',')[[1]][1])}

    if (alt == 'A'){NV=as.integer(strsplit(row[["gt_AU"]], split=',')[[1]][1])}
    if (alt == 'T'){NV=as.integer(strsplit(row[["gt_TU"]], split=',')[[1]][1])}
    if (alt == 'G'){NV=as.integer(strsplit(row[["gt_GU"]], split=',')[[1]][1])}
    if (alt == 'C'){NV=as.integer(strsplit(row[["gt_CU"]], split=',')[[1]][1])}
    if (nchar(ref) > 1 | nchar(alt) > 1 ){NV=as.integer(strsplit(row[["gt_TIR"]], split=',')[[1]][1])}

    ref_alt = paste0(NR, ',', NV)
    ref_alt
}

parse_Strelka = function(vcf, tumour_id, normal_id){

    tb = vcfR::vcfR2tidy(vcf)
    # Strelka tumour_id and normal_id is TUMOR and NORMAL, change to tumour_id/normal_id
    gt_field = tb[["gt"]] %>%
        dplyr::mutate(sample = ifelse(Indiv=="TUMOR", tumour_id, normal_id)) %>%
        dplyr::select(-Indiv)

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

    samples_list = gt_field[["sample"]] %>% unique

    calls = lapply(
        samples_list,
        function(s){
            gt_field_s = gt_field %>% dplyr::filter(sample == s)

            fits = list()
            fits[["sample"]] = s
            fits[["mutations"]] = dplyr::bind_cols(fix_field, gt_field_s)
            mutations = fits[["mutations"]]
            ref_alt = lapply(1:nrow(mutations), function(r){
                retrieve_ref_alt(mutations[r,])
            })

            ref_alt = ref_alt %>% unlist()
            mutations[["ref_alt"]] = ref_alt
            mutations = mutations %>%
                tidyr::separate(ref_alt, into = c('NR', 'NV')) %>%
                dplyr::mutate(NR = as.integer(NR), NV = as.integer(NV)) %>%
                dplyr::mutate(DP = NR+NV) %>%
                dplyr::mutate(VAF = NV/DP) %>%
                dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, everything())

            fits[["mutations"]] = mutations
            fits
    })

    samples = c(normal_id, tumour_id)
    if (normal_id == '[]'){
      samples = c(tumour_id)
    } else{
      samples = c(tumour_id, normal_id)
    }
    names(calls) = samples

    if ("CSQ" %in% tb[['meta']][['ID']]){
    # VEP specific field extraction
    # Take CSQ field names and split by |

        vep_field = tb[['meta']] %>%
            dplyr::filter(ID == "CSQ") %>%
            dplyr::select(Description) %>%
            dplyr::pull()

        tmp_vep_field = strsplit(vep_field, split = "|", fixed = TRUE) %>% unlist()
        vep_field = tmp_vep_field[1:length(tmp_vep_field)-1]

        # Tranform the fix field by splittig the CSQ and select the columns needed
        calls[[tumour_id]][['mutations']] = calls[[tumour_id]][['mutations']] %>%
            dplyr::mutate(CSQ = strsplit(CSQ, ",")) %>%
            tidyr::unnest(CSQ) %>%
            tidyr::separate(CSQ, vep_field, sep = "\\\\|") %>%
            dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything())  #can add other thing, CSQ, HGSP

        if (normal_id != '[]'){
          calls[[normal_id]][['mutations']] = calls[[normal_id]][['mutations']] %>% dplyr::select(-CSQ) %>% dplyr::distinct()
        }
    }
    return(calls)
}

parse_Platypus = function(vcf, tumour_id, normal_id){
    tb = vcfR::vcfR2tidy(vcf)

    gt_field = tb[["gt"]] %>%
        dplyr::mutate(
            DP = as.numeric(gt_NR),
            NV = as.numeric(gt_NV),
            VAF = NV/DP) %>%
        dplyr::rename(sample = Indiv)

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
        dplyr::select(chr, from, to, ref, alt, CSQ, dplyr::everything(), -ChromKey)

    # Extract each sample
    samples_list = gt_field[["sample"]] %>% unique

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
                dplyr::mutate(NR = DP - NV)
            fits
        }
    )

    names(calls) = samples_list
    if (normal_id == '[]'){
      samples = c(tumour_id)
    } else{
      samples = c(tumour_id, normal_id)
    }
    calls = calls[samples]

    if ("CSQ" %in% tb[['meta']][['ID']]){
    # VEP specific field extraction
    # Take CSQ field names and split by |

        vep_field = tb[['meta']] %>%
            dplyr::filter(ID == "CSQ") %>%
            dplyr::select(Description) %>%
            dplyr::pull()

        tmp_vep_field = strsplit(vep_field, split = "|", fixed = TRUE) %>% unlist()
        vep_field = tmp_vep_field[1:length(tmp_vep_field)-1]

        # Tranform the fix field by splittig the CSQ and select the columns needed
        calls[[tumour_id]][['mutations']] = calls[[tumour_id]][['mutations']] %>%
            dplyr::mutate(CSQ = strsplit(CSQ, ",")) %>%
            tidyr::unnest(CSQ) %>%
            tidyr::separate(CSQ, vep_field, sep = "\\\\|") %>%
            dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything())  #can add other thing, CSQ, HGSP

        if (normal_id != '[]'){
          calls[[normal_id]][['mutations']] = calls[[normal_id]][['mutations']] %>% dplyr::select(-CSQ) %>% dplyr::distinct()
        }
    }
    return(calls)
}


parse_TNscope = function(vcf, tumour_id, normal_id){
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
    dplyr::select(chr, from, to, ref, alt, dplyr::everything(),
                  -ChromKey, -ID, -QUAL, -FILTER, -TVAF, -TDP, -NVAF, -NDP, -TAL)

  # Extract sample names
  samples_list = gt_field[["sample"]] %>% unique

  calls = lapply(
    samples_list,
    function(s){
      gt_field_s = gt_field %>% dplyr::filter(sample == s)

      if(nrow(fix_field) != nrow(gt_field_s))
        stop("Mismatch between the VCF fixed fields and the genotypes, will not process this file.")

      fits = list()
      fits[["sample"]] = s
      fits[["mutations"]] = dplyr::bind_cols(fix_field, gt_field_s) %>%
        dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, dplyr::everything())
      fits
    })

  names(calls) = samples_list
  if (normal_id == '[]'){
    samples = c(tumour_id)
  } else{
    samples = c(tumour_id, normal_id)
  }
  calls = calls[samples]

  # check if VCF is annotated with VEP
  if ("CSQ" %in% tb[["meta"]][["ID"]]){
    # VEP specific field extraction
    # Take CSQ field names and split by |

    vep_field = tb[['meta']] %>%
      dplyr::filter(ID == "CSQ") %>%
      dplyr::select(Description) %>%
      dplyr::pull()

    tmp_vep_field = strsplit(vep_field, split = "|", fixed = TRUE) %>% unlist()
    vep_field = tmp_vep_field[1:length(tmp_vep_field)-1]

    # Tranform the fix field by splittig the CSQ and select the columns needed
    calls[[tumour_id]][['mutations']] = calls[[tumour_id]][['mutations']] %>%
      dplyr::mutate(CSQ = strsplit(CSQ, ",")) %>%
      tidyr::unnest(CSQ) %>%
      tidyr::separate(CSQ, vep_field, sep = "\\\\|") %>%
      dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything())  #can add other thing, CSQ, HGSP

    if (normal_id != '[]'){
      calls[[normal_id]][['mutations']] = calls[[normal_id]][['mutations']] %>% dplyr::select(-CSQ) %>% dplyr::distinct()
    }
  }
  return(calls)
}


library(dplyr)
library(tidyr)
library(vcfR)


# Read vcf file
vcf = vcfR::read.vcfR("$vcf")

# Check from which caller the .vcf has been produced
caller_sig <- if (length(vcf@meta)) paste(vcf@meta, collapse = "\n") else ""

if (grepl(pattern = 'TNscope|TNhaplotyper2', x = caller_sig, ignore.case = TRUE)) {
  calls <- parse_TNscope(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample")

} else if (grepl(pattern = 'Mutect', x = caller_sig, ignore.case = TRUE)) {
  calls <- parse_Mutect(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample")

} else if (grepl(pattern = 'strelka', x = caller_sig, ignore.case = TRUE)) {
  calls <- parse_Strelka(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample")

} else if (grepl(pattern = 'Platypus', x = caller_sig, ignore.case = TRUE)) {
  calls <- parse_Platypus(vcf, tumour_id = "$meta.tumour_sample", normal_id = "$meta.normal_sample")

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
