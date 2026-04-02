/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_VIEW } from "..//modules/nf-core/bcftools/view/main"
include { ENSEMBLVEP_VEP } from "../modules/nf-core/ensemblvep/vep/main"
include { FORMATTER as FORMATTER_CNA } from "../subworkflows/local/formatter/main"
include { FORMATTER as FORMATTER_VCF} from "../subworkflows/local/formatter/main"
include { LIFTER } from "../subworkflows/local/lifter/main"
include { ANNOTATE_DRIVER } from "../modules/local/annotate_driver/main"
include { FORMATTER as FORMATTER_RDS} from "../subworkflows/local/formatter/main"
include { QC } from "../subworkflows/local/qc/main"
include { SUBCLONAL_DECONVOLUTION } from "../subworkflows/local/subclonal_deconvolution/main"
include { SIGNATURE_DECONVOLUTION } from "../subworkflows/local/signature_deconvolution/main"

include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TUMOUREVO {

    take:
    input_samplesheet
    fasta
    drivers_table
    vep_cache
    ch_versions

main:
    input = input_samplesheet.map{ meta, vcf, tbi, bam, bai, cna_segs, cna_extra ->
            meta = meta + [id: "${meta.dataset}_${meta.patient}_${meta.tumour_sample}"]
            [meta.dataset + meta.patient, meta, vcf, tbi, bam, bai, cna_segs, cna_extra] }
            .groupTuple()
            .map { id, meta, vcf, tbi, bam, bai, cna_segs, cna_extra ->
                n = vcf.baseName.unique().size()
                [id, meta, vcf, tbi, bam, bai, cna_segs, cna_extra, n ]}
            .transpose()
            .map { id, meta, vcf, tbi, bam, bai, cna_segs, cna_extra, n  ->
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

    if (params.vcf_filter_mutations == true){
        BCFTOOLS_VIEW(input_vcf, [], [], [])
        vcf = BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi)
    } else {
        vcf = input_vcf
    }

    ENSEMBLVEP_VEP(
        vcf,
        params.vep_genome,
        params.vep_species,
        params.vep_cache_version,
        vep_cache,
        fasta,
        ch_extra_files,
    )

    ch_vcf_tbi = ENSEMBLVEP_VEP.out.vcf.join(ENSEMBLVEP_VEP.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    FORMATTER_VCF(ch_vcf_tbi, "vcf")
    FORMATTER_CNA(input_cna, "cna")
    vcf_file = FORMATTER_VCF.out.out_data
    cna_file = FORMATTER_CNA.out.out_data
    ch_versions = ch_versions.mix(FORMATTER_VCF.out.versions)
    ch_versions = ch_versions.mix(FORMATTER_CNA.out.versions)

    join_input = vcf_file.join(input_bam).map{ meta, rds, bam, bai ->
            [ meta, rds, bam, bai ] }
            .branch { meta, rds, bam, bai ->
                to_lift: meta.lifter == true
                multisample: meta.lifter == false
            }

    LIFTER(join_input.to_lift, fasta)
    out_lifter = LIFTER.out.out_data
    ch_versions = ch_versions.mix(LIFTER.out.versions)

    rds_input = join_input.multisample.map{ meta, rds, bam, bai ->
            [meta, rds]
            }

    vcf_rds = rds_input.concat(out_lifter)
    ANNOTATE_DRIVER(vcf_rds.combine(drivers_table))
    ch_versions = ch_versions.mix(ANNOTATE_DRIVER.out.versions)

    in_cnaqc = cna_file.join(ANNOTATE_DRIVER.out.rds)
    QC(in_cnaqc)
    ch_versions = ch_versions.mix(QC.out.versions)

    if (params.filter == true){
        SUBCLONAL_DECONVOLUTION(QC.out.join_cnaqc_PASS)
        SIGNATURE_DECONVOLUTION(QC.out.join_cnaqc_PASS)

    } else {
        SUBCLONAL_DECONVOLUTION(QC.out.join_cnaqc_ALL)
        SIGNATURE_DECONVOLUTION(QC.out.join_cnaqc_ALL)
    }

    ch_versions = ch_versions.mix(SUBCLONAL_DECONVOLUTION.out.versions)
    ch_versions = ch_versions.mix(SIGNATURE_DECONVOLUTION.out.versions)

    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_tumourevo_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

emit:
    versions = ch_collated_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
