//
// FORMATTING SUB-WORKFLOW
//

include { CNAQC2TSV } from '../../../modules/local/cnaqc2tsv/main'
include { CNA2CNAQC } from '../../../modules/local/cna2cnaqc/main'
include { VCF2CNAQC } from '../../../modules/local/vcf2cnaqc/main'


workflow FORMATTER {
    take:
        input
        extension

    main:
        ch_versions = Channel.empty()
        if (extension == "vcf"){
                VCF2CNAQC(input)
                ch_versions = ch_versions.mix(VCF2CNAQC.out.versions)
                out = VCF2CNAQC.out.rds

        } else if (extension == "cna"){
                input_cna = input.map{ meta, cna  ->
                    [ meta, cna[0], cna[1] ]
                }
                CNA2CNAQC(input_cna)
                ch_versions = ch_versions.mix(CNA2CNAQC.out.versions)
                out = CNA2CNAQC.out.rds

        } else if (extension == "rds"){ // for pyclone-vi
                CNAQC2TSV(input)
                ch_versions = ch_versions.mix(CNAQC2TSV.out.versions)
                out = CNAQC2TSV.out.tsv
        }

    emit:
        out_data = out
        versions = ch_versions
}
