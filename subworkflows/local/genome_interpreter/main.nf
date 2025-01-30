//
// GENOME INTERPRETER SUB-WORKFLOW
//

include { XXX } from '../../../modules/local/xxx/main'


workflow GENOME_INTERPRETER {
    take:
        cnaqc
        drivers
        signatures
        clusters
        vcf

    main:
        maf = VCF2MAF(vcf)
        annotated_maf = ADD_GENOMIC_INFO(maf, cnaqc, drivers,clusters, signatures) // a process that attaches information from copy numbers, drivers, subclonal and signature deconvolution to the maf files
        maf_plot = MAFPLOT(maf)
        xxx

    emit:
        xxx
}
