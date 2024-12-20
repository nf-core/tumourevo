//
// DRIVER_ANNOTATION SUB-WORKFLOW
// 

//include { BUILD_REFERENCE } from '../../../modules/local/build_reference/main.nf'
//include { DNDSCV } from '../../../modules/local/dndscv/main.nf'
include { ANNOTATE_DRIVER } from '../../../modules/local/annotate_driver/main.nf'
//include { JOIN_ANNOTATION } from '../../../modules/local/join_annotation/main.nf'


workflow DRIVER_ANNOTATION {
    take:        
        rds
        driver_list
        //cds
        //genome
    
    main:
        //if (params.dndscv_refcds_rda) {
        //    rda = Channel.from(file(params.dndscv_refcds_rda, checkIfExists: true))
        //    in_dnds = rds.map{meta, data ->
        //        meta = meta + [id: "${meta.dataset}"]
        //        sample = meta.tumour_sample
        //        patient = meta.patient
        //        [meta.subMap('dataset', 'id'), data, patient, sample]}
        //        | groupTuple
        //    dndscv_ch = in_dnds.combine(rda)

        //} else {
        //    // create reference for variant annotation, meta here is not needed
        //    BUILD_REFERENCE(cds, genome)
        //    rda = BUILD_REFERENCE.out.dnds_reference
        //    dndscv_ch = rds.combine(driver_list).combine(rda)
        //}
        
        ANNOTATE_DRIVER(rds.combine(driver_list))
        //DNDSCV(dndscv_ch) // add dndsCV statistics and columns "known_driver" (based on driver list, default IntoGen) and "potential_driver"

        //DNDSCV.out.annot_rds.view()
        //JOIN_ANNOTATION(ANNOTATE_DRIVER.out.rds.combine(DNDSCV.out.annot_rds))

    emit:
        //dnds_annot = DNDSCV.out.annot_rds
        //dnds_res = DNDSCV.out.dnds_rds
        annot_rds = ANNOTATE_DRIVER.out.rds

        //full_annot = JOIN_ANNOTATION.out.rds
}
