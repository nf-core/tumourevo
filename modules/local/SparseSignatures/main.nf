process SPARSE_SIGNATURES {
  tag "$meta.id"
  container = 'docker://lvaleriani/sparsesignature:v2'
  

  input:

    tuple val(meta), path(tsv_join,  stageAs: '*.tsv')

  output:
 
    tuple val(meta), path("*_cv_means_mse.rds"), emit: signatures_cv_rds
    tuple val(meta), path("*_best_params_config.rds"), emit: signatures_bestConf_rds
    tuple val(meta), path("*_nmf_Lasso_out.rds"), emit: signatures_nmfOut_rds
    tuple val(meta), path("*_plot_all.rds"), emit: signatures_plot_rds
    tuple val(meta), path("*_plot_all.pdf"), emit: signatures_plot_pdf
                            
  script:

    def args                         = task.ext.args 
    def prefix                       = task.ext.prefix                                    ?: "${meta.id}"
    def K                            = args!='' && args.K                                 ? "$args.K" : ""
    def background_signature         = args!='' && args.background_signature              ? "$args.background_signature" : ""
    def beta                         = args!='' && args.beta                              ? "$args.beta" : ""
    def normalize_counts             = args!='' && args.normalize_counts                  ? "$args.normalize_counts" : ""
    def nmf_runs                     = args!='' && args.nmf_runs                          ? "$args.nmf_runs" : ""
    def iterations                   = args!='' && args.iterations                        ? "$args.iterations" : ""
    def max_iterations_lasso         = args!='' && args.max_iterations_lasso              ? "$args.max_iterations_lasso" : ""
    def num_processes                = args!='' && args.num_processes                     ? "$args.num_processes" : ""
    def cross_validation_entries     = args!='' && args.cross_validation_entries          ? "$args.cross_validation_entries" : ""
    def cross_validation_repetitions = args!='' && args.cross_validation_repetitions      ? "$args.cross_validation_repetitions" : ""
    def cross_validation_iterations  = args!='' && args.cross_validation_iterations       ? "$args.cross_validation_iterations" : ""
    def lambda_values_alpha          = args!='' && args.lambda_values_alpha               ? "$args.lambda_values_alpha" : ""
    def lambda_values_beta           = args!='' && args.lambda_values_beta                ? "$args.lambda_values_beta" : ""
    def lambda_rate_alpha            = args!='' && args.lambda_rate_alpha                 ? "$args.lambda_rate_alpha" : ""
    def verbose                      = args!='' && args.verbose                           ? "$args.verbose" : ""
    def seed                         = args!='' && args.seed                              ? "$args.seed" : ""

    """
    #!/usr/bin/env Rscript

    library(SparseSignatures)
    library(ggplot2)
    library(stringr)
    library(patchwork)
    library(dplyr)

    source("$moduleDir/getters.R")

    n_procs = parse(text="$num_processes")
    if (n_procs == "all"){
	n_procs = eval("Inf")
    }

    patients_tsv = strsplit("$tsv_join", " ")[[1]]
    tables = lapply(patients_tsv, FUN = function(p_table){
        read.delim(p_table, sep = "\\t", header=T) %>% 
          mutate(across(everything(), as.character)) 
      } 
    ) 
    multisample_table = dplyr::bind_rows(tables)


    #Extract input data information
    input_data <- multisample_table[,c("Indiv","chr","from","to","ref","alt")]
    input_data <- setNames(input_data, c("sample","chrom","start","end","ref","alt"))
    input_data[["end"]] <- input_data[["start"]]
    input_data <- input_data %>% mutate(start = as.integer(start), end = as.integer(end))

    #Generate the patient vs mutation count matrix from mutation data
    #Load reference human-genome specification.
    #The user must select, among the available choices, the reference genome consistent with the mutation dataset.

    load_genome <- function(genome, input_data) {
      if (genome == "GRCh37") {
        library(BSgenome.Hsapiens.1000genomes.hs37d5)
        bsg <- BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5
        input_data[["chrom"]] <- substring(input_data[["chrom"]], 4, 5)

      } else if (genome == "GRCh38") {
        library(BSgenome.Hsapiens.UCSC.hg38)
        bsg <- BSgenome.Hsapiens.UCSC.hg38

        # Leave 'chrom' unchanged for GRCh38
      }
      return(list(bsg = bsg, input_data = input_data))
    }

    data_list = load_genome("$params.genome", input_data)
    bsg <- data_list[["bsg"]]
    input_data <- data_list[["input_data"]] 

    
    mut_counts = SparseSignatures::import.trinucleotides.counts(data=input_data, reference=bsg)

    #Load a reference SBS5 background signature from COSMIC
    data(background)

    #Estimate the initial values of beta
    starting_betas = SparseSignatures::startingBetaEstimation(x = mut_counts, 
                                                              K = eval(parse(text="$K")),
                                                              background_signature  = background)


    #Find the optimal number of signatures and sparsity level: rely on cross-validation
    # higher number of CV repetitions corresponds to more accurate parameter estimates

    cv_out = SparseSignatures::nmfLassoCV(
      x = mut_counts,
      K = eval(parse(text="$K")),  
      starting_beta = starting_betas,
      background_signature = background, 
      normalize_counts = as.logical("$normalize_counts"), 
      nmf_runs = as.integer("$nmf_runs"), 
      lambda_values_alpha = eval(parse(text="$lambda_values_alpha")), 
      lambda_values_beta = eval(parse(text="$lambda_values_beta")),
      cross_validation_entries = as.numeric("$cross_validation_entries"), 
      cross_validation_iterations = as.integer("$cross_validation_iterations"), 
      cross_validation_repetitions = as.integer("$cross_validation_repetitions"), 
      iterations = as.integer("$iterations"), 
      max_iterations_lasso = as.integer("$max_iterations_lasso"), 
      num_processes = n_procs, 
      verbose = as.logical("$verbose"),
      seed = as.integer("$seed")
    ) 


    #Analyze the mean squared error results averaging over cross-validation repetitions
    cv_mses <- cv_out[["grid_search_mse"]][1,,]
    cv_means_mse <- matrix(sapply(cv_mses, FUN = mean),
                          nrow = dim(cv_mses)[1]
    ) 

    dimnames(cv_means_mse) <- dimnames(cv_mses)

    #Find the combination of parameters that yields the lowest MSE
    min_ii <- which(cv_means_mse == min(cv_means_mse, na.rm = TRUE), arr.ind = TRUE)
    min_Lambda_beta <- rownames(cv_means_mse)[min_ii[1]]
    min_Lambda_beta <- as.numeric(gsub("_Lambda_Beta", "", min_Lambda_beta))
    min_K <- colnames(cv_means_mse)[min_ii[2]] 
    min_K <- as.numeric(gsub("_Signatures", "", min_K))
    best_params_config <- data.frame(min_K, min_Lambda_beta)

    saveRDS(object = cv_means_mse, file = paste0("$prefix", "_cv_means_mse.rds"))
    saveRDS(object = best_params_config, file = paste0("$prefix", "_best_params_config.rds")) 

    #Discovering the signatures within the dataset: NMF Lasso
    #Compute the signatures for the best configuration.

    nmf_Lasso_out = SparseSignatures::nmfLasso(
      x = mut_counts,
      K = min_K,
      beta = eval(parse(text="$beta")), 
      background_signature = background, 
      normalize_counts = as.logical("$normalize_counts"),
      lambda_rate_alpha = eval(parse(text="$lambda_rate_alpha")), 
      lambda_rate_beta = min_Lambda_beta,
      iterations = as.integer("$iterations"), 
      max_iterations_lasso = as.integer("$max_iterations_lasso"), 
      verbose = as.logical("$verbose")
    ) 

    saveRDS(object = nmf_Lasso_out, file =  paste0("$prefix", "_nmf_Lasso_out.rds"))

    #Signature visualization
    signatures = nmf_Lasso_out\$beta
    plot_signatures <- SparseSignatures::signatures.plot(beta=signatures, xlabels=FALSE)

    plot_exposure = nmf_Lasso_out\$alpha %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var="PatientID") %>% 
      tidyr::pivot_longer(cols=!"PatientID", names_to="Signatures", values_to="Exposures") %>% 
    
      ggplot() +
      geom_bar(aes(x=PatientID, y=Exposures, fill=Signatures), 
              position="stack", stat="identity") +
      theme(axis.text.x=element_text(angle=90,hjust=1),
            panel.background=element_blank(),
            axis.line=element_line(colour="black"))

    plt_all = patchwork::wrap_plots(plot_exposure, plot_signatures, ncol=2) + patchwork::plot_annotation(title = "$meta.id")
    ggplot2::ggsave(plot = plt_all, filename = paste0("$prefix", "_plot_all.pdf"), width = 210, height = 297, units="mm", dpi = 200)
    saveRDS(object = plt_all, file = paste0("$prefix", "_plot_all.rds"))

    """ 
 }
