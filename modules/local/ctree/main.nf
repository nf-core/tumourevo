process CTREE {
    tag "$meta.id"
    // container='file:///fast/cdslab/ebusca00/singularity/cdslab.sif'
    container = 'docker://elenabuscaroli/ctree:latest'

    input:

        tuple val(meta), path(ctree_input)

    output:
        tuple val(meta), path("**ctree_{VIBER,MOBSTERh,pyclonevi}.rds"), emit: ctree_rds, optional: true
        tuple val(meta), path("**ctree_{VIBER,MOBSTERh,pyclonevi}_plots.rds"), emit: ctree_plots_rds, optional: true
        tuple val(meta), path("**REPORT_plots_ctree_{VIBER,MOBSTERh,pyclonevi}.rds"), emit: ctree_report_rds, optional: true
        tuple val(meta), path("**REPORT_plots_ctree_{VIBER,MOBSTERh,pyclonevi}.pdf"), emit: ctree_report_pdf, optional: true
        tuple val(meta), path("**REPORT_plots_ctree_{VIBER,MOBSTERh,pyclonevi}.png"), emit: ctree_report_png, optional: true

    script:

    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sspace_cutoff = args!="" && args.sspace_cutoff ? "$args.sspace_cutoff" : ""
    def n_sampling = args!="" && args.n_sampling ? "$args.n_sampling" : ""
    def store_max = args!="" && args.store_max ? "$args.store_max" : ""

    """
    #!/usr/bin/env Rscript

    library(ctree)
    library(dplyr)
    library(VIBER)
    library(mobster)
    library(ggplot2)

    outdir = ""

    add_dummy_driver = function(input_table, variant_colname, is_driver_colname) {
        add_driver = FALSE
        idx = dplyr::case_when(
            "cluster" %in% colnames(input_table) ~ which(input_table[["cluster"]] != "Tail")[1],
            .default = 1
        )
        if (!variant_colname %in% colnames(input_table) | !is_driver_colname %in% colnames(input_table)) {
            input_table = input_table %>% dplyr::mutate(!!is_driver_colname:=FALSE, !!variant_colname:=NA)
            add_driver = TRUE
        } else if (all(input_table[[is_driver_colname]]==FALSE)) {
            # idx = which(input_table[["cluster"]] != "Tail")[1]
            add_driver = TRUE
        }
        if (add_driver) {
            input_table[idx, is_driver_colname] = TRUE
            input_table[idx, variant_colname] = ""
        }
        return(input_table)
    }

    initialize_ctree_obj_pyclone = function(ctree_input) {
    # ctree_input = add_dummy_driver(ctree_input, variant_colname="variantID", is_driver_colname="is.driver")

        driver_cluster = unique(ctree_input[which(ctree_input["is.driver"]==TRUE),c("cluster")])
        # the CCF table must report CCF values for each cluster and sample
        # cluster | nMuts | is.driver | is.clonal | sample1 | sample2 | ...
        CCF_table = ctree_input %>%
            dplyr::select(sample_id, cluster, nMuts, is.driver, is.clonal, CCF) %>%

        dplyr::mutate(is.driver=replace(is.driver, is.driver=="", "FALSE")) %>%
        dplyr::mutate(is.driver=as.logical(is.driver)) %>%

        dplyr::filter(cluster!="Tail") %>%
        dplyr::mutate(cluster=as.character(cluster)) %>%

        dplyr::group_by(cluster) %>%
        dplyr::mutate(is.driver=any(is.driver)) %>%
        dplyr::filter(any(CCF>0)) %>%
        dplyr::ungroup() %>% unique() %>%
        tidyr::pivot_wider(names_from="sample_id", values_from="CCF", values_fill=0)

        # the driver table must contain patient and variant IDs and report clonality and driver status
        # patientID | variantID | is.driver | is.clonal | cluster | sample1 | sample2 | ...
        drivers_table = ctree_input %>%
            dplyr::filter(cluster %in% CCF_table[["cluster"]]) %>%
            dplyr::mutate(is.driver=as.logical(is.driver)) %>%
            dplyr::mutate(cluster=as.character(cluster)) %>%
            dplyr::select(patientID, sample_id, variantID, cluster, is.driver, is.clonal, CCF) %>%
            dplyr::filter(is.driver==TRUE) %>%
            dplyr::mutate(variantID=replace(variantID, is.na(variantID), "")) %>%
            tidyr::pivot_wider(names_from="sample_id", values_from="CCF", values_fill=0)

        samples = unique(ctree_input[["sample_id"]])  # if multisample, this is a list
        patient = unique(ctree_input[["patientID"]])

        CCF_table = add_dummy_driver(CCF_table, variant_colname="variantID", is_driver_colname="is.driver")


        if (nrow(drivers_table)==0) {
            drivers_table = CCF_table %>%
            dplyr::filter(is.driver) %>%
            dplyr::select(-nMuts) %>%
            dplyr::mutate(patientID=patient)
        }

        ctree_init = list("CCF_table"=CCF_table,
                            "drivers_table"=drivers_table,
                            "samples"=samples,
                            "patient"=patient)
        return(ctree_init)
    }

    if ( grepl(".rds\$", tolower("$ctree_input")) ) {
        best_fit = readRDS("$ctree_input")
        do_fit = TRUE

        if (!"data" %in% names(best_fit)) { best_fit[["data"]] =  data.frame(row.names=1:best_fit[["N"]])}

        ## viber
        if (class(best_fit) == "vb_bmm") {
            fn_name = VIBER::get_clone_trees
            subclonal_tool = "VIBER"

        # If only one sample inference ctree is not working with VIBER
        if (ncol(best_fit[["x"]]) <= 2) {
            cli::cli_alert_warning("Object of class {class(best_fit)} with only one sample is not supported by the function `get_clone_trees`!")
            do_fit = FALSE
        }

        best_fit[["data"]] = add_dummy_driver(best_fit[["data"]], variant_colname="gene", is_driver_colname="driver")
        }

        ## mobster
        if (class(best_fit) == "dbpmm") {
            fn_name = mobster::get_clone_trees
            subclonal_tool = "MOBSTERh"
            sample_id = unique(best_fit[["data"]][["sample_id"]])
            outdir = paste0(sample_id, "/")

        best_fit[["data"]] = add_dummy_driver(best_fit[["data"]], variant_colname="driver_label", is_driver_colname="is_driver")
        }

        if (class(best_fit) %in% c("vb_bmm", "dbpmm") & do_fit) {
            # VIBER or MOBSTER object
            trees = fn_name(x = best_fit,
                            sspace.cutoff = as.integer("$sspace_cutoff"),
                            n.sampling = as.integer("$n_sampling"),
                            store.max = as.integer("$store_max"))
        } else {
            cli::cli_alert_warning("Object of class {class(best_fit)} not supported.")
            do_fit = FALSE
        }

    } else {
        do_fit = TRUE
        subclonal_tool = "pyclonevi"
        input_table = read.csv("$ctree_input", sep="\t")
        data_ctree = initialize_ctree_obj_pyclone(input_table)
        trees = ctrees(CCF_clusters = data_ctree[["CCF_table"]],
                        drivers = data_ctree[["drivers_table"]],
                        samples = data_ctree[["samples"]],
                        patient = data_ctree[["patient"]],
                        sspace.cutoff = as.integer("$sspace_cutoff"),
                        n.sampling = as.integer("$n_sampling"),
                        store.max = as.integer("$store_max"))
    }

    if (do_fit & !is.null(trees)) {
        ctree_output = paste0("ctree_", subclonal_tool)
        if (!dir.exists(outdir)) dir.create(outdir, recursive=T)

        # plot the best tree
        top_phylo = plot(trees[[1]])

        # save rds and plots
        saveRDS(object=trees, file=paste0(outdir, "$prefix", "_", ctree_output, ".rds"))
        saveRDS(object=top_phylo, file=paste0(outdir, "$prefix", "_", ctree_output, "_plots.rds"))

        # Save report plot
        phylos = ggplot2::ggplot()
        if (length(trees) > 1) phylos = lapply(trees[2:min(length(trees), 3)], plot)
        phylos = ggpubr::ggarrange(plotlist=phylos)

        ccf = ctree::plot_CCF_clusters(trees[[1]])
        info_transfer = ctree::plot_information_transfer(trees[[1]])
        clone_size = ctree::plot_clone_size(trees[[1]])

        report_fig = ggpubr::ggarrange(plotlist=list(ccf, info_transfer, top_phylo, clone_size, phylos), nrow=3, ncol=2)

        saveRDS(object=report_fig, file=paste0(outdir, "$prefix", "_REPORT_plots_", ctree_output, ".rds"))
        ggplot2::ggsave(plot=report_fig, filename=paste0(outdir, "$prefix", "_REPORT_plots_", ctree_output, ".pdf"), height=297, width=210, units="mm", dpi=200)
        ggplot2::ggsave(plot=report_fig, filename=paste0(outdir, "$prefix", "_REPORT_plots_", ctree_output, ".png"), height=297, width=210, units="mm", dpi=200)
    }

    """
}
