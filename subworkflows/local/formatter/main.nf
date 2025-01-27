//
// FORMATTING SUB-WORKFLOW
//

include { RDS_PROCESSING } from '../../../modules/local/CNAqc2tsv/main'
include { CNA_PROCESSING } from '../../../modules/local/cna2CNAqc/main'
include { VCF_PROCESSING } from '../../../modules/local/vcf2CNAqc/main'


workflow FORMATTER {
    take:
        input
        extension

    main:

        if (extension == "vcf"){
                VCF_PROCESSING(input)
                out = VCF_PROCESSING.out.rds

        } else if (extension == "cna"){
                input_cna = input.map{ meta, cna  ->
                    [ meta, cna[0], cna[1] ]
                }
                CNA_PROCESSING(input_cna)
                out = CNA_PROCESSING.out.rds

        } else if (extension == "rds"){ // for pyclone-vi
                RDS_PROCESSING(input)
                out = RDS_PROCESSING.out.tsv
        }

    emit:
        out
}

