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


# Auxiliary functions #####

get_sample = function(m_cnaqc_obj,
                    sample,
                    which_obj) {
    if (class(m_cnaqc_obj) != "m_cnaqc") {
        cli::cli_abort(c(
        "mCNAqc object was build without keeping original CNAqc objects",
        "x" = "It is not possible to retrieve the required samples"
        ))
    }

    consented_obj = c("shared", "original")

    if ((which_obj %in% consented_obj) == FALSE) {
        cli::cli_abort("{.var which_obj} must be one of 'shared' or 'original'")
    }
    # define the element names
    if (which_obj == "original") {
        type = "original_cnaqc_objc"
        # check if the original cnaqc obj exist
        check_or = any(names(m_cnaqc_obj) == type)
        if (check_or == FALSE) {
        cli::cli_abort(c(
            "mCNAqc object was build without keeping original CNAqc objects",
            "x" = "It is not possible to retrieve the required samples"
        ))
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
        cli::cli_abort(c(
        "Must provide a {.field m_cnaqc} object",
        "x" = "{.var x} is a {.cls {class(x)}}"
        ))
    }
}

get_mCNAqc_stats = function(m_cnaqc_obj) {
    stats = m_cnaqc_obj[["m_cnaqc_stats"]]
    return(stats)
}


# Script #####

library(CNAqc)
library(mobster)
library(dplyr)
library(ggplot2)

patientID = description = "$meta.patient"
samples = strsplit(x = "$meta.tumour_sample", ",") %>% unlist()  # list of samples

## read mCNAqc object
if ( grepl(".rds\$", tolower("$rds_join")) ) {
    obj = readRDS("$rds_join")
    if (class(obj) == "m_cnaqc") {
        original = obj %>% get_sample(sample=samples, which_obj="original")
        input_table = lapply(names(original),
                                function(sample_name)
                                original[[sample_name]] %>%
                                    # keep only mutations on the diploid karyotype
                                    CNAqc::subset_by_segment_karyotype("1:1") %>%
                                    CNAqc::Mutations() %>%
                                    dplyr::mutate(sample_id=sample_name)
                                ) %>% dplyr::bind_rows()
    } else {
        cli::cli_abort("Object of class {class($rds_join)} not supported.")
    }
} else {
    input_table = read.csv("$rds_join")
}

## Function to run a single mobster fit
run_mobster_fit = function(joint_table, descr) {
    # get input table for the single patient
    inp_tb = joint_table %>%
        dplyr::filter(VAF < 1) %>%
        # dplyr::mutate(VAF=replace(VAF, VAF==0, 1e-7)) %>%
        dplyr::filter(VAF!=0) %>%
        dplyr::filter(karyotype=="1:1")
        # dplyr::rename(variantID=gene) %>%
        # dplyr::rename(is.driver=is_driver) %>%
        # dplyr::rename(tumour_content=purity) %>%
        # dplyr::rename(is_driver=is.driver)
        # dplyr::rename(is_driver=is.driver, driver_label=variantID)

# mobster_fit(x = inp_tb,
#             K = eval(parse(text="$K")),
#             samples = as.integer("$samples"),
#             init = "$init",
#             tail = eval(parse(text="$tail")),
#             epsilon = as.numeric("$epsilon"),
#             maxIter = as.integer("$maxIter"),
#             fit.type = "$fit_type",
#             seed = as.integer("$seed"),
#             model.selection = "$model_selection",
#             trace = as.logical("$trace"),
#             parallel = as.logical("$parallel"),
#             pi_cutoff = as.numeric("$pi_cutoff"),
#             N_cutoff = as.integer("$n_cutoff"),
#             auto_setup = eval(parse(text="$auto_setup")),
#             silent = as.logical("$silent"),
#             description = descr)
# }

mobster_fit(x = inp_tb,
            K = eval(parse(text=opt[["K"]])),
            samples = as.integer(opt[["samples"]]),
            init = opt[["init"]],
            tail = eval(parse(text=opt[["tail"]])),
            epsilon = as.numeric(opt[["epsilon"]]),
            maxIter = as.integer(opt[["maxIter"]]),
            fit.type = opt[["fit_type"]],
            seed = as.integer(opt[["seed"]]),
            model.selection = opt[["model_selection"]],
            trace = as.logical(opt[["trace"]]),
            parallel = as.logical(opt[["parallel"]]),
            pi_cutoff = as.numeric(opt[["pi_cutoff"]]),
            N_cutoff = as.integer(opt[["n_cutoff"]]),
            auto_setup = eval(parse(text=opt[["auto_setup"]])),
            silent = as.logical(opt[["silent"]]),
            description = descr)
}

lapply(samples, function(sample_name) {

    fit = run_mobster_fit(joint_table=input_table %>% dplyr::filter(sample_id == !!sample_name),
                        descr=description)

    best_fit = fit[["best"]]
    plot_fit = plot(best_fit)

    saveRDS(object=fit, file=paste0(opt[["prefix"]], "_mobsterh_st_fit.rds"))
    saveRDS(object=best_fit, file=paste0(opt[["prefix"]], "_mobsterh_st_best_fit.rds"))
    saveRDS(object=plot_fit, file=paste0(opt[["prefix"]], "_mobsterh_st_best_fit_plots.rds"))

    # save report plots
    report_fig = mobster::plot_model_selection(fit)
    saveRDS(report_fig, file=paste0(opt[["prefix"]], "_REPORT_plots_mobster.rds"))
    ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_REPORT_plots_mobster.pdf"), height=210, width=210, units="mm", dpi = 200)
    ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_REPORT_plots_mobster.png"), height=210, width=210, units="mm", dpi = 200)
})


# version export
f = file("versions.yml","w")
cnaqc_version = sessionInfo()\$otherPkgs\$CNAqc\$Version
mobster_version = sessionInfo()\$otherPkgs\$mobster\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    mobster:", mobster_version), f)
close(f)
