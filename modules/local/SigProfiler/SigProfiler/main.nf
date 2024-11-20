process SIGPROFILER {
    tag "$meta.id"
    container = 'docker://katiad/sigprofiler:latest'
  
    input:
       tuple val(meta), path(tsv_list, stageAs: '*.tsv')
       path(genome_path)

    output:
       tuple val(meta), path("results/*"), emit: sigprofiler_results
      
    script:
    
      def args                              = task.ext.args                                 ?: ''
      def prefix                            = task.ext.prefix                               ?: "${meta.id}"
      def reference_genome                  = args!='' && args.reference_genome             ? "$args.reference_genome" : ""
      def exome                             = args!='' && args.exome                        ? "$args.exome" : ""
      def input_type                        = args!='' && args.input_type                   ? "$args.input_type" : ""
      def context_type                      = args!='' && args.context_type                 ? "$args.context_type" : ""
      def minimum_signatures                = args!='' && args.minimum_signatures           ? "$args.minimum_signatures" : ""
      def maximum_signatures                = args!='' && args.maximum_signatures           ? "$args.maximum_signatures" : ""
      def nmf_replicates                    = args!='' && args.nmf_replicates               ? "$args.nmf_replicates" : ""
      def resample                          = args!='' && args.resample                     ? "$args.resample" : ""
      def seeds                             = args!='' && args.seeds                        ? "$args.seeds" : ""
      def matrix_normalization              = args!='' && args.matrix_normalization         ? "$args.matrix_normalization" : ""
      def nmf_init                          = args!='' && args.nmf_init                     ? "$args.nmf_init" : "random"
      def min_nmf_iterations                = args!='' && args.min_nmf_iterations           ? "$args.min_nmf_iterations" : ""
      def max_nmf_iterations                = args!='' && args.max_nmf_iterations           ? "$args.max_nmf_iterations" : ""
      def nmf_test_conv                     = args!='' && args.nmf_test_conv                ? "$args.nmf_test_conv" : ""
      def nmf_tolerance                     = args!='' && args.nmf_tolerance                ? "$args.nmf_tolerance" : ""
      def cpu                               = args!='' && args.cpu                          ? "$args.cpu" : ""
      def stability                         = args!='' && args.stability                    ? "$args.stability" : ""
      def min_stability                     = args!='' && args.min_stability                ? "$args.min_stability" : ""
      def combined_stability                = args!='' && args.combined_stability           ? "$args.combined_stability" : ""  
      def cosmic_version                    = args!='' && args.cosmic_version               ? "$args.cosmic_version" : ""
      def make_decomposition_plots          = args!='' && args.make_decomposition_plots     ? "$args.make_decomposition_plots" : ""
      def collapse_to_SBS96                 = args!='' && args.collapse_to_SBS96            ? "$args.collapse_to_SBS96" : ""
      def get_all_signature_matrices        = args!='' && args.get_all_signature_matrices   ? "$args.get_all_signature_matrices" : ""
      def export_probabilities              = args!='' && args.export_probabilities         ? "$args.export_probabilities": ""   
      def bed_file                          = args!='' && args.bed_file                     ? "$args.bed_file" : ""   
      def chrom_based                       = args!='' && args.chrom_based                  ? "$args.chrom_based" : ""          
      def plot                              = args!='' && args.plot                         ? "$args.plot" : ""       
      def tsb_stat                          = args!='' && args.tsb_stat                     ? "$args.tsb_stat" : ""                 
      def seqInfo                           = args!='' && args.seqInfo                      ? "$args.seqInfo" : ""             
      def cushion                           = args!='' && args.cushion                      ? "$args.cushion" : ""             
      def precision                         = args!='' && args.precision                    ? "$args.precision" : ""                  
      def gpu                               = args!='' && args.gpu                          ? "$args.gpu" : ""      
      def batch_size                        = args!='' && args.batch_size                   ? "$args.batch_size" : ""                   
      def allow_stability_drop              = args!='' && args.allow_stability_drop         ? "$args.allow_stability_drop" : ""                                         
  
    
      """
      #!/usr/bin/env python3
     
      import os
      import shutil
      import pandas as pd
      import multiprocessing
      from SigProfilerExtractor import sigpro as sig
      from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
      #from utils_sigprofiler import input_processing, process_tsv_join      

      if __name__ == '__main__':
          
          dataset_id = "$meta.dataset"
          input_path = os.path.join(dataset_id)

          if not os.path.exists(input_path):
              os.mkdir(input_path)
      
          output_path = os.path.join("output", "SBS", f"{dataset_id}_SBS96_all")
         
          
          # input data preprocessing

          def process_tsv_join(tsv_list):
             patients_tsv = tsv_list.split()
             # Read each file into a pandas DataFrame and ensure all columns are of type 'string'
             tables = []
             for p_table in patients_tsv:
                 df = pd.read_csv(p_table, sep='\\t', dtype=str)
                 tables.append(df)
             multisample_table = pd.concat(tables, ignore_index=True)
             return multisample_table

          def input_processing(data, dataset_id, genome = '$referece_genome'):
             new_columns = {'Project': dataset_id, 'Genome': genome, 'Type': "SOMATIC", 'mut_type': "SNP"}
             df = data.assign(**new_columns)
             df['chr'] = df['chr'].astype(str).str[3:]
             df = df.rename(columns={'Indiv': 'Sample', 'chr': 'chrom', 'from': 'pos_start', 'to': 'pos_end'})
             df["ID"] = df["Sample"]
             df = df.loc[:, ['Project', 'Sample', 'ID', 'Genome', 'mut_type', 'chrom', 'pos_start', 'pos_end', 'ref', 'alt', 'Type']]
             return df
    
          input_tsv_join = process_tsv_join("$tsv_list")

          input_data = input_processing(input_tsv_join, dataset_id, '$reference_genome')

          # saving input matrix to txt
          input_data.to_csv(f"{input_path}/input_data.txt", sep="\\t", index=False, header=True)

          # mutation's counts matrix generation
          input_matrix = matGen.SigProfilerMatrixGeneratorFunc(
                  project = dataset_id, 
                  reference_genome = '$reference_genome', 
                  path_to_input_files = input_path,
                  exome = bool("$exome"),
                  bed_file = "$bed_file", 
                  chrom_based = bool("$chrom_based"), 
                  plot = bool("$plot"), 
                  tsb_stat = bool("$tsb_stat"), 
                  seqInfo = bool("$seqInfo"), 
                  cushion = int("$cushion"),
                  volume = "./")

          full_input_data_path = os.path.join(input_path, output_path)

          # Perform model fitting
          sig.sigProfilerExtractor(input_type = "$input_type", 
                                   output = "results", 
                                   input_data = full_input_data_path,  
                                   reference_genome = '$reference_genome', 
                                   opportunity_genome = '$reference_genome',
                                   context_type = "$context_type",  
                                   exome = bool("$exome"),
                                   minimum_signatures = int("$minimum_signatures"),  
                                   maximum_signatures = int("$maximum_signatures"), 
                                   nmf_replicates = int("$nmf_replicates"),
                                   resample = bool("$resample"),
                                   seeds= "$seeds",
                                   matrix_normalization = "$matrix_normalization", 
                                   nmf_init = "$nmf_init", 
                                   precision = "$precision",
                                   min_nmf_iterations = int("$min_nmf_iterations"), 
                                   max_nmf_iterations = int("$max_nmf_iterations"),
                                   nmf_test_conv = int("$nmf_test_conv"), 
                                   nmf_tolerance = float("$nmf_tolerance"), 
                                   cpu = int("$cpu"),
                                   gpu = bool("$gpu"),
                                   batch_size = int("$batch_size"),
                                   stability = float("$stability"),
                                   min_stability = float("$min_stability"),
                                   combined_stability = float("$combined_stability"),
                                   allow_stability_drop = bool("$allow_stability_drop"),
                                   cosmic_version = float("$cosmic_version"),
                                   make_decomposition_plots = bool("$make_decomposition_plots"), 
                                   collapse_to_SBS96 = bool("$collapse_to_SBS96"), 
                                   get_all_signature_matrices = bool("$get_all_signature_matrices"),
                                   export_probabilities = bool("$export_probabilities"))

          # save the output results
          dest_dir = "$prefix/"
          source_dir = "results/"
          shutil.copytree(source_dir, dest_dir, dirs_exist_ok=True)
   """
}
