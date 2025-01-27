process MOBSTERh {
    tag "$meta.id"
    label "process_single"
    container = 'docker://elenabuscaroli/mobster:version1.0'

    input:
        tuple val(meta), path(rds_join) // rds from JOIN_CNAQC

    output:
        tuple val(meta), path("*_mobsterh_st_fit.rds"), emit: mobster_rds
        tuple val(meta), path("*_mobsterh_st_best_fit.rds"), emit: mobster_best_rds
        tuple val(meta), path("*_plots.rds"), emit: mobster_plots_rds
        tuple val(meta), path("*_REPORT_plots_mobster.rds"), emit: mobster_report_rds
        tuple val(meta), path("*_REPORT_plots_mobster.pdf"), emit: mobster_report_pdf
        tuple val(meta), path("*_REPORT_plots_mobster.png"), emit: mobster_report_png
        path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def K = args!="" && args.K ? "$args.K" : ""
    def samples = args!="" && args.samples ? "$args.samples" : ""
    def init = args!="" && args.init ? "$args.init" : ""
    def tail = args!="" && args.tail ? "$args.tail" : ""
    def epsilon = args!="" && args.epsilon ? "$args.epsilon" : ""
    def maxIter = args!="" && args.maxIter ? "$args.maxIter" : ""
    def fit_type = args!="" && args.fit_type ? "$args.fit_type" : ""
    def seed = args!="" && args.seed ? "$args.seed" : ""
    def model_selection = args!="" && args.model_selection ? "$args.model_selection" : ""
    def trace = args!="" && args.trace ? "$args.trace" : ""
    def parallel = args!="" && args.parallel ? "$args.parallel" : ""
    def pi_cutoff = args!="" && args.pi_cutoff ? "$args.pi_cutoff" : ""
    def n_cutoff = args!="" && args.n_cutoff ? "$args.n_cutoff" : ""
    def auto_setup = args!="" && args.auto_setup ? "$args.auto_setup" : ""
    def silent = args!="" && args.silent ? "$args.silent" : ""

    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
    library(CNAqc)
    library(mobster)
    library(dplyr)
    library(ggplot2)
    source("$moduleDir/getters.R")
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

    mobster_fit(x = inp_tb,
                K = eval(parse(text="$K")),
                samples = as.integer("$samples"),
                init = "$init",
                tail = eval(parse(text="$tail")),
                epsilon = as.numeric("$epsilon"),
                maxIter = as.integer("$maxIter"),
                fit.type = "$fit_type",
                seed = as.integer("$seed"),
                model.selection = "$model_selection",
                trace = as.logical("$trace"),
                parallel = as.logical("$parallel"),
                pi_cutoff = as.numeric("$pi_cutoff"),
                N_cutoff = as.integer("$n_cutoff"),
                auto_setup = eval(parse(text="$auto_setup")),
                silent = as.logical("$silent"),
                description = descr)
    }

    lapply(samples, function(sample_name) {
    #  dir.create(outDir_sample, recursive = TRUE)

        fit = run_mobster_fit(joint_table=input_table %>% dplyr::filter(sample_id == !!sample_name),
                            descr=description)

        best_fit = fit[["best"]]
        plot_fit = plot(best_fit)

        saveRDS(object=fit, file=paste0("$prefix", "_mobsterh_st_fit.rds"))
        saveRDS(object=best_fit, file=paste0("$prefix", "_mobsterh_st_best_fit.rds"))
        saveRDS(object=plot_fit, file=paste0("$prefix", "_mobsterh_st_best_fit_plots.rds"))

        # save report plots
        report_fig = mobster::plot_model_selection(fit)
        saveRDS(report_fig, file=paste0("$prefix", "_REPORT_plots_mobster.rds"))
        ggplot2::ggsave(plot=report_fig, filename=paste0("$prefix", "_REPORT_plots_mobster.pdf"), height=210, width=210, units="mm", dpi = 200)
        ggplot2::ggsave(plot=report_fig, filename=paste0("$prefix", "_REPORT_plots_mobster.png"), height=210, width=210, units="mm", dpi = 200)
    })


    # version export
    f = file("versions.yml","w")
    cnaqc_version = sessionInfo()\$otherPkgs\$CNAqc\$Version
    mobster_version = sessionInfo()\$otherPkgs\$mobster\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    CNAqc:", cnaqc_version), f)
    writeLines(paste("    mobster:", mobster_version), f)
    close(f)
    """

}
