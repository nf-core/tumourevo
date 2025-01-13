process VIBER {
    tag "$meta.id"
    container = 'docker://elenabuscaroli/viber:version0.1'

    input:
    tuple val(meta), path(rds_join), val(tumour_samples) //rds from either JOIN_CNAQC or JOIN_FIT, should be always grouped

    output:
    tuple val(meta), path("*_viber_best_st_fit.rds"), emit: viber_rds
    tuple val(meta), path("*_viber_best_st_heuristic_fit.rds"), emit: viber_heuristic_rds
    tuple val(meta), path("*_${plot1}"), emit: viber_plots_rds
    tuple val(meta), path("*_${plot2}"), emit: viber_heuristic_plots_rds
    tuple val(meta), path("*_REPORT_plots_viber.rds"), emit: viber_report_rds
    tuple val(meta), path("*_REPORT_plots_viber.pdf"), emit: viber_report_pdf
    tuple val(meta), path("*_REPORT_plots_viber.png"), emit: viber_report_png

script:
    // viber fit params
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?:"${meta.id}"
    def K = args!="" && args.K ? "$args.K" : ""
    def alpha_0 = args!="" && args.alpha_0 ? "$args.alpha_0" : ""
    def a_0 = args!="" && args.a_0 ? "$args.a_0" : ""
    def b_0 = args!="" && args.b_0 ? "$args.b_0" : ""
    def maxIter = args!="" && args.maxIter ? "$args.maxIter" : ""
    def epsilon_conv = args!="" && args.epsilon_conv ? "$args.epsilon_conv" : ""
    def samples = args!="" && args.samples ? "$args.samples" : ""
    def q_init = args!="" && args.q_init ? "$args.q_init" : ""
    def trace = args!="" && args.trace ? "$args.trace" : ""
    // viber filter params
    def binomial_cutoff = args!="" && args.binomial_cutoff ? "$args.binomial_cutoff" : ""
    def dimensions_cutoff = args!="" && args.dimensions_cutoff ? "$args.dimensions_cutoff" : ""
    def pi_cutoff = args!="" && args.pi_cutoff ? "$args.pi_cutoff" : ""
    def re_assign = args!="" && args.re_assign ? "$args.re_assign" : ""
    // def mode = args.mode ?  "$args.mode" : ""
    def n_samples = tumour_samples.size()

    if (n_samples==1) {
        plot1 = "viber_best_st_mixing_plots.rds"
        plot2 = "viber_best_st_heuristic_mixing_plots.rds"
    } else {
        plot1 = "viber_best_st_fit_plots.rds"
        plot2 = "viber_best_st_heuristic_fit_plots.rds"
    }

    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

    library(VIBER)
    library(dplyr)
    library(tidyverse)
    library(ggplot2)
    library(CNAqc)
    source("$moduleDir/getters.R")

    patientID = "$meta.patient"
    samples = substr("$tumour_samples", 2, nchar("$tumour_samples")-1)
    samples = strsplit(samples, ", ")[[1]]
    print("$tumour_samples")
    print("$rds_join")
    print(samples)

    if ( grepl(".rds\$", tolower("$rds_join")) ) {
        input_obj = readRDS("$rds_join")
        if (class(input_obj) == "m_cnaqc") {
            shared = input_obj %>% get_sample(sample=samples, which_obj="shared")
            joint_table = lapply(names(shared),
                            function(sample_name)
                            shared[[sample_name]] %>%
                                # keep only mutations on the diploid karyotype
                                CNAqc::subset_by_segment_karyotype("1:1") %>%
                                CNAqc::Mutations() %>%
                                dplyr::mutate(sample_id=sample_name)
                            ) %>% dplyr::bind_rows()
    } else {
        cli::cli_alert_warning("Object of class {class(input_obj)} not supported.")
        return()
        }
    }

    print("Subset joint done")

    ## TODO : add drivers to `input_tab`

    ## Read input joint table
    input_tab = joint_table %>%
    dplyr::mutate(VAF=replace(VAF, VAF==0, 1e-7))

    ## Convert the input table into longer format
    reads_data = input_tab %>%
    dplyr::select(chr, from, ref, alt, NV, DP, VAF, sample_id) %>%
    tidyr::pivot_wider(names_from="sample_id",
                        values_from=c("NV","DP","VAF"), names_sep=".")

    ## Extract DP (depth)
    dp = reads_data %>%
    # dplyr::filter(mutation_id %in% non_tail) %>% ## this step should be managed before by other module
    dplyr::select(dplyr::starts_with("DP")) %>%
    dplyr::mutate(dplyr::across(.cols=dplyr::everything(),
                                .fns=function(x) replace(x, is.na(x), 0))) %>%
    dplyr::rename_all(function(x) stringr::str_remove_all(x,"DP."))

    ## Extract NV (alt_counts)
    nv = reads_data %>%
    # dplyr::filter(mutation_id %in% non_tail) %>% ## this step should be managed before by other module
    dplyr::select(dplyr::starts_with("NV")) %>%
    dplyr::mutate(dplyr::across(.cols=dplyr::everything(),
                                .fns=function(x) replace(x, is.na(x), 0))) %>%
    dplyr::rename_all(function(x) stringr::str_remove_all(x,"NV."))

    # Standard fit
    viber_K = eval(parse(text="$K"))
    viber_K[which.min(viber_K)] = 2
    st_fit = VIBER::variational_fit(nv, dp,
                                    K=viber_K,
                                    data=reads_data,
                                    # %>% dplyr::filter(mutation_id %in% non_tail)
                                    alpha_0=as.numeric("$alpha_0"),
                                    a_0=as.integer("$a_0"),
                                    b_0=as.integer("$b_0"),
                                    max_iter=as.integer("$maxIter"),
                                    epsilon_conv=as.numeric("$epsilon_conv"),
                                    samples=as.integer("$samples"),
                                    q_init="$q_init",
                                    trace=as.logical("$trace"),
                                    description=""
                                    )

    best_fit = best_fit_heuristic = st_fit

    # If all clusters are removed -> keep the origianl best fit
    tryCatch(expr = {
    # Apply the heuristic
    best_fit_heuristic = VIBER::choose_clusters(st_fit,
        binomial_cutoff=as.numeric("$binomial_cutoff"),
        dimensions_cutoff=as.integer("$dimensions_cutoff"),
        pi_cutoff=as.numeric("$pi_cutoff"),
        re_assign=as.logical("$re_assign")
        )
    }, error = function(e) {
    print(e)
    best_fit_heuristic <<- st_fit
    } )

    # Save fits
    saveRDS(best_fit, file=paste0("$prefix", "_viber_best_st_fit.rds"))
    saveRDS(best_fit_heuristic, file = paste0("$prefix", "_viber_best_st_heuristic_fit.rds"))

    # Save plots
    if ("$n_samples" >1) { #mutlisample mode on
    print("multisample mode on")
    plot_fit = plot(best_fit)
    plot_fit_heuristic = plot(best_fit_heuristic)

    saveRDS(plot_fit, file=paste0("$prefix", "_viber_best_st_fit_plots.rds"))
    saveRDS(plot_fit_heuristic, file=paste0("$prefix", "_viber_best_st_heuristic_fit_plots.rds"))
    } else {
    plot_fit_mixing = plot_mixing_proportions(best_fit)
    plot_fit_mixing_heuristic = plot_mixing_proportions(best_fit_heuristic)

    saveRDS(plot_fit_mixing, file=paste0("$prefix", "_viber_best_st_mixing_plots.rds"))
    saveRDS(plot_fit_mixing_heuristic, file=paste0("$prefix", "_viber_best_st_heuristic_mixing_plots.rds"))
    }

    # Save report plot
    n_samples = ncol(best_fit[["x"]]) - 1
    marginals = multivariate = ggplot()

    try(expr = {marginals <<- VIBER::plot_1D(best_fit)} )
    #try(expr = {multivariate = plot(best_fit) %>% patchwork::wrap_plots()} )
    #top_p = patchwork::wrap_plots(marginals, multivariate, design=ifelse(n_samples>2, "ABB", "AAB"))

    try(expr = {multivariate = plot(best_fit)})
    try(expr = {multivariate = ggpubr::ggarrange(plotlist = multivariate)})
    top_p = ggpubr::ggarrange(plotlist = list(marginals, multivariate), widths=ifelse(n_samples>2, c(1,2), c(2,1)))

    mix_p = VIBER::plot_mixing_proportions(best_fit)
    binom = VIBER::plot_peaks(best_fit)
    elbo = VIBER::plot_ELBO(best_fit)
    #bottom_p = patchwork::wrap_plots(mix_p, binom, elbo, design="ABBBC")
    bottom_p = ggpubr::ggarrange(plotlist = list(mix_p, binom, elbo), widths = c(1,2,1), nrow = 1)

    #report_fig = patchwork::wrap_plots(top_p, bottom_p, design=ifelse(n_samples>2, "AAAB", "AAB"))
    report_fig = ggpubr::ggarrange(top_p, bottom_p, nrow = 2, heights = ifelse(n_samples>2, c(3,1), c(2,1)))
    saveRDS(report_fig, file=paste0("$prefix", "_REPORT_plots_viber.rds"))
    ggplot2::ggsave(plot=report_fig, filename=paste0("$prefix", "_REPORT_plots_viber.pdf"), height=210, width=210, units="mm", dpi = 200)
    ggplot2::ggsave(plot=report_fig, filename=paste0("$prefix", "_REPORT_plots_viber.png"), height=210, width=210, units="mm", dpi = 200)

    """

    stub:
    """
    echo "${task.process}:" > versions.yml
    echo ' CNAqc: 1.0.0' >> versions.yml
    echo ' viber: 0.1.0' >> versions.yml
    """
}
