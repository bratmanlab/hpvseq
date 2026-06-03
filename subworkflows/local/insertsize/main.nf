//
// Coverage for quantification on properly-paired reads
//
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PROPERLYPAIRED     } from '../../../modules/nf-core/samtools/view'
//include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index'
include { PICARD_COLLECTINSERTSIZEMETRICS  } from '../../../modules/nf-core/picard/collectinsertsizemetrics'

workflow INSERTSIZE {
    take:
    input           // channel: [ val(meta), path(bam), path(bai) ]

    main:

    //
    // Extract properly-paired reads (-f 2)
    //
    SAMTOOLS_VIEW_PROPERLYPAIRED (
        input,
        [[:], [], []],
        [[:], []],
        [[:], []],
        ""
    )
    ch_properlypaired = SAMTOOLS_VIEW_PROPERLYPAIRED.out.bam
//    SAMTOOLS_INDEX(ch_properlypaired)
//    ch_insertsize = ch_properlypaired.join(SAMTOOLS_INDEX.out.index)

    PICARD_COLLECTINSERTSIZEMETRICS (
        ch_properlypaired
    )

    emit:
    metrics      = PICARD_COLLECTINSERTSIZEMETRICS.out.metrics 
}

