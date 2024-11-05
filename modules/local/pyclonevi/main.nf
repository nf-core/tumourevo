process PYCLONEVI {
    tag "$meta.id"
    container = 'https://depot.galaxyproject.org/singularity/pyclone-vi%3A0.1.3--pyhca03a8a_0'

    input:
      tuple val(meta), path(rds_join), val(tumour_samples)

    output:
      tuple val(meta), path("*_cluster_table.csv"), emit: ctree_input
      tuple val(meta), path("*.tsv"), emit: pyclone_input
      tuple val(meta), path("*_all_fits.h5"), emit: pyclone_all_fits
      tuple val(meta), path("*_best_fit.txt"), emit: pyclone_best_fit

    script:
      def args = task.ext.args ?: ''
      // def prefix = task.ext.prefix ?:"${meta.id}_remove_tail_$args.remove_tail" 
      def prefix = task.ext.prefix ?:"${meta.id}" 
      def n_cluster_arg = args.n_cluster ? "$args.n_cluster" : ""
      def density_arg = args.density ? "$args.density" : ""
      def n_grid_point_arg = args.n_grid_point ? "$args.n_grid_point" : ""
      def n_restarts_arg = args.n_restarts ? "$args.n_restarts" : ""
      sampleID_string = tumour_samples.join(" ")

      """
      # format the input table in order to be pyclone compliant
      python3 $moduleDir/pyclone_utils.py create_pyclone_input $rds_join $meta.patient ${prefix}_pyclone_input_all_samples.tsv

      colnames=\$(head -n 1 ${prefix}_pyclone_input_all_samples.tsv)

      column_number=\$(echo -e "\$colnames" | awk -v col_name="sample_id" 'BEGIN { FS = "\t" } {
          for (i = 1; i <= NF; i++) {
              if (\$i == col_name) {
                  print i
                  exit
              }
          }
          exit 1  # Exit with error if column not found
      }')

      echo -e "\$colnames" > ${prefix}_pyclone_input.tsv
      for i in $sampleID_string;
        do awk '\$'\$column_number' == "'"\$i"'"' ${prefix}_pyclone_input_all_samples.tsv >> ${prefix}_pyclone_input.tsv;
      done

      pyclone-vi fit -i ${prefix}_pyclone_input.tsv -o ${prefix}_all_fits.h5 -c $n_cluster_arg -d $density_arg --num-grid-points $n_grid_point_arg --num-restarts $n_restarts_arg
      pyclone-vi write-results-file -i ${prefix}_all_fits.h5 -o ${prefix}_best_fit.txt

      python3 $moduleDir/pyclone_ctree.py --joint ${prefix}_pyclone_input.tsv --best_fit ${prefix}_best_fit.txt --ctree_input ${prefix}_cluster_table.csv
      """
}
