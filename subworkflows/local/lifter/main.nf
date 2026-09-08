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
        def ch_out      = channel.empty()
        def ch_versions = channel.empty()

        def bam = data.map{ meta, _rds, _bam, _bai ->
            [meta, _bam]
        }

        def rds = data.map{ meta, _rds, _bam, _bai ->
            [meta, _rds]
        }

        def all_rds = data.map{ meta, _rds, _bam, _bai ->
            def new_meta = meta + [id: "${meta.dataset}_${meta.patient}"]
            [new_meta.subMap('dataset', 'patient', 'id', 'normal_sample'), _rds]
        }.groupTuple()

        GET_POSITIONS_ALL(all_rds)
        ch_versions = ch_versions.mix(GET_POSITIONS_ALL.out.versions)

        def all_pos = GET_POSITIONS_ALL.out.all_pos.transpose().map{ meta, _rds ->
            [_rds]
        }

        GET_POSITIONS_REL(rds.combine(all_pos))
        ch_versions = ch_versions.mix(GET_POSITIONS_REL.out.versions)

        def in_pileup = bam.join(GET_POSITIONS_REL.out.bed, by: [0])

        BCFTOOLS_MPILEUP(in_pileup, fasta, false)
        ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)

        def ch_join = rds.join(BCFTOOLS_MPILEUP.out.vcf, by: [0])

        JOIN_POSITIONS(ch_join.combine(all_pos))
        ch_versions = ch_versions.mix(JOIN_POSITIONS.out.versions)
        ch_out = JOIN_POSITIONS.out.rds

    emit:
        out_data = ch_out
        versions = ch_versions

}
