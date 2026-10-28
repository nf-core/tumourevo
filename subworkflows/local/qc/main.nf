//
// QC SUB-WORKFLOW
//

include { TINC } from '../../../modules/nf-core/tinc/main'
include { CNAQC as CNAQC } from '../../../modules/nf-core/cnaqc/main'
include { JOIN_CNAQC } from '../../../modules/local/join_cnaqc/main'
include { CNAQC as CNAQC_BY_CHR } from '../../../modules/nf-core/cnaqc/main'
include { JOIN_CNAQC as JOIN_CNAQC_BY_CHR } from '../../../modules/local/join_cnaqc/main'

workflow QC {
    take:
        input

    main:
        ch_versions = Channel.empty()

        plot_rds_tinc = null
        rds_tinc = null
        pdf_tinc = null
        csv_tinc = null

        plot_cnaqc_by_chr_qc = null
        plot_cnaqc_by_chr_qc_pdf = null
        rds_cnaqc_by_chr = null

        if (params.tools && params.tools.split(',').contains('tinc')) {
            TINC(input)
            ch_versions = ch_versions.mix(TINC.out.versions_tinc)
            plot_rds_tinc = TINC.out.plot_rds
            rds_tinc = TINC.out.rds
            pdf_tinc = TINC.out.plot_pdf
            csv_tinc = TINC.out.tinc_csv
        }

        // set the input for cnaqc
        input_cnaqc = input.map{meta, cna, snv ->
                def sample =  meta.tumour_sample
                [meta, snv, cna, sample]
        }

        // check if the flag is set to true or not

        if(params.qc_chr == true) {

            CNAQC_BY_CHR(input_cnaqc)

            plot_cnaqc_by_chr_qc_pdf = CNAQC_BY_CHR.out.plot_pdf_qc_by_chr
            plot_cnaqc_by_chr_qc = CNAQC_BY_CHR.out.plot_qc_by_chr
            rds_cnaqc_by_chr = CNAQC_BY_CHR.out.qc_rds
            ch_versions = ch_versions.mix(CNAQC_BY_CHR.out.versions)

            // emit results from whole genome analysis
            rds_cnaqc = CNAQC_BY_CHR.out.qc_rds
            plot_cnaqc_data_rds = CNAQC_BY_CHR.out.data_plot_rds
            plot_cnaqc_qc_rds = CNAQC_BY_CHR.out.qc_plot_rds
            plot_cnaqc_data = CNAQC_BY_CHR.out.plot_pdf_data
            plot_cnaqc_qc = CNAQC_BY_CHR.out.plot_pdf_qc

            // join by chromosomes
            in_join_cnaqc = CNAQC_BY_CHR.out.qc_by_chr_rds.map{ meta, rds ->
                def sample = meta.tumour_sample
                meta = meta + [id: "${meta.dataset}_${meta.patient}"]
                [meta.subMap('dataset', 'patient', 'id'), rds, sample]}
                .groupTuple()

            JOIN_CNAQC_BY_CHR(in_join_cnaqc)
            join_cnaqc_ALL = JOIN_CNAQC_BY_CHR.out.rds_all
            join_cnaqc_PASS = JOIN_CNAQC_BY_CHR.out.rds_pass

            ch_versions = ch_versions.mix(JOIN_CNAQC_BY_CHR.out.versions)

        } else {

            CNAQC(input_cnaqc)
            ch_versions = ch_versions.mix(CNAQC.out.versions)

            rds_cnaqc = CNAQC.out.qc_rds
            plot_cnaqc_data_rds = CNAQC.out.data_plot_rds
            plot_cnaqc_qc_rds = CNAQC.out.qc_plot_rds
            plot_cnaqc_data = CNAQC.out.plot_pdf_data
            plot_cnaqc_qc = CNAQC.out.plot_pdf_qc

            in_join_cnaqc = CNAQC.out.qc_rds.map{ meta, rds ->
                def sample = meta.tumour_sample
                meta = meta + [id: "${meta.dataset}_${meta.patient}"]
                [meta.subMap('dataset', 'patient', 'id'), rds, sample]}
                .groupTuple()

            JOIN_CNAQC(in_join_cnaqc)
            join_cnaqc_ALL = JOIN_CNAQC.out.rds_all
            join_cnaqc_PASS = JOIN_CNAQC.out.rds_pass

            ch_versions = ch_versions.mix(JOIN_CNAQC.out.versions)
        }

    emit:
        rds_cnaqc
        plot_cnaqc_data_rds
        plot_cnaqc_qc_rds
        plot_cnaqc_data
        plot_cnaqc_qc

        plot_rds_tinc
        rds_tinc
        pdf_tinc
        csv_tinc

        plot_cnaqc_by_chr_qc
        plot_cnaqc_by_chr_qc_pdf
        rds_cnaqc_by_chr

        join_cnaqc_ALL
        join_cnaqc_PASS

        versions = ch_versions
}
