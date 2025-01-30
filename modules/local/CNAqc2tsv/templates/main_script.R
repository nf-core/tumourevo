#!/usr/bin/env Rscript

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)

# Auxiliary functions #####

get_sample = function(m_cnaqc_obj, sample, which_obj) {
    if (class(m_cnaqc_obj) != "m_cnaqc") {
        wrong_class_all = class(m_cnaqc_obj)

        cli::cli_abort(
            c("cnaqc_objs must be a {.field m_cnaqc} object",
            "x" = "{.var m_cnaqc_obj} is a {.cls {class(m_cnaqc_obj)}}")
        )
    }

    consented_obj = c("shared", "original")
    if ((which_obj %in% consented_obj) == FALSE) {
        cli::cli_abort("{.var which_obj} must be one of {.val shared} or {.val original}")
    }

    # define the element names
    if (which_obj == "original") {
        type = "original_cnaqc_objc"
        # check if the original cnaqc obj exist
        check_or = any(names(m_cnaqc_obj) == type)
        if (check_or == FALSE) {
            cli::cli_abort(c("mCNAqc object was build without keeping original CNAqc objects"),
                        "x" = "It is not possible to retrieve the required samples")
        } else {
            cli::cli_h1("Retrieving original {.cls CNAqc} objects")
            cnaqc_samples = m_cnaqc_obj[[type]][sample]
        }

    } else {
        type = "cnaqc_obj_new_segmentation"
        cli::cli_h1("Retrieving {.cls CNAqc} objects with the new segmentation")
        cnaqc_samples = m_cnaqc_obj[[type]][sample]
    }
    return(cnaqc_samples)
}

get_sample_name = function(x) {
    if (class(x) == "m_cnaqc") {
        lapply(x[["cnaqc_obj_new_segmentation"]], function(y) {
            y[["sample"]]
        }) %>% unlist() %>% unname()
    } else if (class(x) == "cnaqc") {
        x[["sample"]]
    } else {
        wrong_class_all = class(x)
        cli::cli_abort(
            c("must provide a {.field m_cnaqc} object",
            "x" = "{.var x} is a {.cls {class(x)}}")
        )
    }
}

get_mCNAqc_stats = function(m_cnaqc_obj){
    stats = m_cnaqc_obj[["m_cnaqc_stats"]]
    return(stats)
}

# Script ####

library(dplyr)
library(CNAqc)

multi_cnaqc = readRDS(file = "$rds_join")
mutations_multisample = get_sample(m_cnaqc_obj = multi_cnaqc,sample = get_sample_name(multi_cnaqc),
                                which_obj = "original")
multisample_jointTable = list()

for (s in get_sample_name(multi_cnaqc)){
    purity = mutations_multisample[[s]][["purity"]]
    multisample_jointTable[[s]] = mutations_multisample[[s]][["mutations"]] %>%
        dplyr::mutate(purity=purity)
}

joint_table = bind_rows(multisample_jointTable)

write.table(joint_table, file = paste0(opt[["prefix"]], "_joint_table.tsv"), append = F, quote = F, sep = "\t", row.names = FALSE)

# version export
f = file("versions.yml","w")
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
cnaqc_version = sessionInfo()\$otherPkgs\$CNAqc\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
close(f)

