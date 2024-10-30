//
// QC SUB-WORKFLOW
//

include { TINC } from '../../../modules/local/tinc/main'
include { CNAQC } from '../../../modules/local/CNAqc/main'
include { JOIN_CNAQC } from '../../../modules/local/join_CNAqc/main'


workflow QC {
    take: 
        input

    main:
        TINC(input)
        contamination = TINC.out.tinc_csv
            .splitCsv( header: true )
            .map{ meta, csv -> 
            meta = meta + [id: "${meta.dataset}_${meta.patient}"] 
            normal_contamination = csv.normal_contamination_flag
            [meta.subMap('dataset', 'patient', 'id'), normal_contamination ]}
            .unique()
            | groupTuple
            | map{ meta, normal_contamination -> 
                normal_contamination = normal_contamination.max()
                [meta, normal_contamination]
            }
        
        // contamination.view()
        
        // TINC_output = QC.out.rds_join.map{  meta, rds, sample -> 
        //     [meta, rds, sample] }
        //     .branch { meta, rds, sample -> 
        //             pass: meta.normal_contamination == '0'
        //             not_pass: meta.normal_contamination == '1'
        // }
        
        // TINC_output = pass_qc.pass
        
        CNAQC(input)

        in_cnaqc = CNAQC.out.qc_rds.map{ meta, rds -> 
            sample = meta.tumour_sample
            meta = meta + [id: "${meta.dataset}_${meta.patient}"]
            [meta, rds, sample]}
            | groupTuple
            | join(contamination)
        
        pass_tinc = in_cnaqc.map{  meta, rds, sample -> 
                [meta, rds, sample] }
                .branch { meta, rds, sample -> 
                        pass: meta.normal_contamination == 0
                        not_pass: meta.normal_contamination == 1
                }
        
        pass_tinc.pass.view()
        // in_cnaqc.view()
        //.join(TINC.out.tinc_csv)

        // in_cnaqc 

        // in_join = in_cnaqc.map{ meta, rds, csv -> 
        //     meta = meta + [id: "${meta.dataset}_${meta.patient}"]
        //     sample = meta.tumour_sample
        //     [meta.subMap('dataset', 'patient', 'id'), rds, csv, sample]}
        //     | groupTuple
        
        // join_cnaqc_out = JOIN_CNAQC(in_join)
        // rds_join = join_cnaqc_out.join(contamination).map{
        //     meta, rds, sample, normal_contamination -> 
        //         meta = meta + [normal_contamination: "${normal_contamination}"]
        //         [meta, rds, sample]
        // }
        // join_cnaqc_out

    
    emit:
        rds_cnaqc = CNAQC.out.qc_rds
        plot_cnaqc_rds = CNAQC.out.plot_rds
        plot_cnaqc_data = CNAQC.out.plot_pdf_data
        plot_cnaqc_qc = CNAQC.out.plot_pdf_qc

        // plot_rds_tinc = TINC.out.plot_rds
        // rds_tinc = TINC.out.rds
        // pdf_tinc = TINC.out.plot_pdf
        // csv_tinc = TINC.out.tinc_csv

        // join_cnaqc_out
        // rds_join
}
