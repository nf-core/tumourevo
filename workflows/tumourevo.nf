include { VCF_ANNOTATE_ENSEMBLVEP } from '../subworkflows/nf-core/vcf_annotate_ensemblvep/main'
include { FORMATTER as FORMATTER_CNA } from "${baseDir}/subworkflows/local/formatter/main"
include { FORMATTER as FORMATTER_VCF} from "${baseDir}/subworkflows/local/formatter/main"
include { LIFTER } from "${baseDir}/subworkflows/local/lifter/main"
include { DRIVER_ANNOTATION } from "${baseDir}/subworkflows/local/annotate_driver/main"
include { FORMATTER as FORMATTER_RDS} from "${baseDir}/subworkflows/local/formatter/main"
include { QC } from "${baseDir}/subworkflows/local/QC/main"
include { SUBCLONAL_DECONVOLUTION } from "${baseDir}/subworkflows/local/subclonal_deconvolution/main"
include { SIGNATURE_DECONVOLUTION } from "${baseDir}/subworkflows/local/signature_deconvolution/main"

workflow TUMOUREVO {

  take:
    input_samplesheet
    fasta
    drivers_table
    vep_cache
  
  main:
    input = input_samplesheet.map{ meta, vcf, tbi, bam, bai, cna_segs, cna_extra -> 
            meta = meta + [id: "${meta.dataset}_${meta.patient}_${meta.tumour_sample}"]
            [meta.dataset + meta.patient, meta, vcf, tbi, bam, bai, cna_segs, cna_extra] }
            | groupTuple 
            | map { id, meta, vcf, tbi, bam, bai, cna_segs, cna_extra -> 
                n = vcf.baseName.unique().size()
                [id, meta, vcf, tbi, bam, bai, cna_segs, cna_extra, n ]}
            | transpose
            | map { id, meta, vcf, tbi, bam, bai, cna_segs, cna_extra, n  ->
                if (n > 1 && bam){
                    meta = meta + [lifter:true]
                } else {
                    meta = meta + [lifter:false]
                }
                [meta, vcf, tbi, bam, bai, cna_segs, cna_extra]
            }


    input_vcf = input.map{ meta, vcf, tbi, bam, bai, cna_segs, cna_extra  -> 
            [ meta, vcf, tbi ]
            }

    input_cna = input.map{ meta, vcf, tbi, bam, bai, cna_segs, cna_extra  -> 
            [ meta, [cna_segs, cna_extra] ]
            }

    input_bam = input.map{ meta, vcf, tbi, bam, bai, cna_segs, cna_extra  -> 
            [ meta, bam, bai ]
            } 

    ch_extra_files = []
    VCF_ANNOTATE_ENSEMBLVEP(input_vcf, 
                            fasta,
                            params.vep_genome,
                            params.vep_species,
                            params.vep_cache_version,
                            vep_cache,
                            ch_extra_files)
    vcf_file = FORMATTER_VCF(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi, "vcf")
    FORMATTER_CNA(input_cna, "cna")

    join_input = vcf_file.join(input_bam).map{ meta, rds, bam, bai -> 
            [ meta, rds, bam, bai ] }
            .branch { meta, rds, bam, bai ->
                to_lift: meta.lifter == true
                multisample: meta.lifter == false
            }

    out_lifter = LIFTER(join_input.to_lift, fasta) 
    rds_input = join_input.multisample.map{ meta, rds, bam, bai -> 
            [meta, rds]
            }
    vcf_rds = rds_input.concat(out_lifter)

    //cds = []
    annotation = DRIVER_ANNOTATION(vcf_rds, drivers_table) //, cds, fasta)
    
    //annotation = vcf_rds
    cna_out = FORMATTER_CNA.out

    in_cnaqc = cna_out.join(annotation)
    QC(in_cnaqc)

    // pass_qc = QC.out.rds_join.map{  meta, rds, sample -> 
    //             [meta, rds, sample] }
    //             .branch { meta, rds, sample -> 
    //                     pass: meta.normal_contamination == '0'
    //                     not_pass: meta.normal_contamination == '1'
    //             }
    
    SUBCLONAL_DECONVOLUTION(QC.out.join_cnaqc_out)
    SIGNATURE_DECONVOLUTION(QC.out.join_cnaqc_out)
}
