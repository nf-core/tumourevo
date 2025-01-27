//
// LIFTER SUB-WORKFLOW
//

include { GET_POSITIONS_ALL } from "../../../modules/local/get_positions/main"
include { GET_POSITIONS_REL } from "../../../modules/local/get_positions/main_rel"
include { BCFTOOLS_MPILEUP } from '../../../modules/nf-core/bcftools/mpileup/main'
include { JOIN_POSITIONS } from "../../../modules/local/join_positions/main"


workflow LIFTER {
    take:
        data
        fasta

    main:
        out = Channel.empty()

        bam = data.map{ meta, rds, bam, bai ->
            [meta, bam]
        }

        rds = data.map{ meta, rds, bam, bai ->
                [meta, rds]
        }

        all_rds = data.map{ meta, rds, bam, bai ->
            meta = meta + [id: "${meta.dataset}_${meta.patient}"]
            [meta.subMap('dataset', 'patient', 'id', 'normal_sample'), rds]}
            | groupTuple

        GET_POSITIONS_ALL(all_rds)
        all_pos = GET_POSITIONS_ALL.out.all_pos.transpose().map{ meta, rds ->
                [rds]
        }

        GET_POSITIONS_REL(rds.combine(all_pos))
        in_pileup = bam.join(GET_POSITIONS_REL.out.bed, by: [0])

        BCFTOOLS_MPILEUP(in_pileup, fasta, false)
        join = rds.join(BCFTOOLS_MPILEUP.out.vcf, by:[0])

        JOIN_POSITIONS(join.combine(all_pos))
        out = JOIN_POSITIONS.out.rds

    emit:
        out

}
