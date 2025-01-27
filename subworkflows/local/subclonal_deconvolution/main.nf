// SUBCLONAL DECONVOLUTION WORKFLOW

include { MOBSTERh } from "../../../modules/local/mobsterh/main"
include { VIBER } from "../../../modules/local/viber/main"
include { PYCLONEVI } from "../../../modules/local/pyclonevi/main"
include { FORMATTER } from "../../../subworkflows/local/formatter/main"
include { CTREE as CTREE_MOBSTERh } from "../../../modules/local/ctree/main"
include { CTREE as CTREE_PYCLONEVI } from "../../../modules/local/ctree/main"
include { CTREE as CTREE_VIBER } from "../../../modules/local/ctree/main"

workflow SUBCLONAL_DECONVOLUTION {
    take:
    rds_join // tuple val(meta), path("*.rds"), val(tumour_samples), emit: rds

    main:
    mobster_pdf = null
    ctree_mobster_pdf = null
    viber_pdf = null
    ctree_viber_pdf = null
    pyclone_fits = null
    pyclone_best = null
    pyclone_table = null
    ctree_pyclone_pdf = null

    if (params.tools && params.tools.split(",").contains("mobster")) {
        joinCNAqc = rds_join.transpose().map{ meta, rds, sample ->
            meta = meta + ["tumour_sample": sample, "id":"${meta.dataset}_${meta.patient}_${sample}"]
            [meta, rds]}
        MOBSTERh(joinCNAqc)
        // in_join = MOBSTERh.out.mobster_best_rds.map{ meta, rds ->
        //     meta = meta + [id: "${meta.dataset}_${meta.patient}"]
        //     sample = meta.tumour_sample
        //     [meta.subMap("dataset", "patient", "id"), rds, sample]}
        //     | groupTuple
        // input_joint_fit = rds_join.map{meta, rds, sample-> [meta, rds]}
        // input_joint_fit.join(in_join).view()
        // rds_join = JOINT_FIT(input_joint_fit.join(in_join))

        CTREE_MOBSTERh(MOBSTERh.out.mobster_best_rds)

        mobster_pdf = MOBSTERh.out.mobster_report_pdf
        ctree_mobster_pdf = CTREE_MOBSTERh.out.ctree_report_pdf

    } // else if (params.remove_tail && !params.remove_tail.contains("never")){
    //     error "None method for tail deconvolution specified"
    // }

    if (params.tools && params.tools.split(",").contains("viber")) {
        VIBER(rds_join)
        CTREE_VIBER(VIBER.out.viber_rds)

        viber_pdf = VIBER.out.viber_report_pdf
        ctree_viber_pdf = CTREE_VIBER.out.ctree_report_pdf
    }

    if (params.tools && params.tools.split(",").contains("pyclone-vi")) {
        FORMATTER(rds_join, "rds")
        PYCLONEVI(FORMATTER.out.out)
        CTREE_PYCLONEVI(PYCLONEVI.out.ctree_input)
        pyclone_fits = PYCLONEVI.out.pyclone_all_fits
        pyclone_best = PYCLONEVI.out.pyclone_best_fit
        pyclone_table = FORMATTER.out.out
        ctree_pyclone_pdf = CTREE_PYCLONEVI.out.ctree_report_pdf
    }

    emit:
    pyclone_fits
    pyclone_best
    ctree_pyclone_pdf
    viber_pdf
    ctree_viber_pdf
    mobster_pdf
    ctree_mobster_pdf
    pyclone_table
}
