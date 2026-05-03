//
// CHECKMATE on SMASH using DEDUP BAM 
//

include { PICARD_MARKDUPLICATES as DEDUP  } from '../../../modules/nf-core/picard/markduplicates'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index'
include { SMASH                   } from '../../../modules/local/smash'

workflow CHECKMATE {
    take:
    bam         // channel: [ val(meta), path(bam) ]
    fasta       // channel: [ val(meta2), path(fasta), path(fai) ]

    main:

    //
    // DE-DUP 
    //
    DEDUP (
        bam,
        fasta
    )
    ch_dedup = DEDUP.out.bam
    SAMTOOLS_INDEX(ch_dedup)
    ch_bam = ch_dedup.join(SAMTOOLS_INDEX.out.index)
    ch_smash = ch_bam
    .map { meta, bam, bai -> [ bam, bai ]}
    .collect()
    .map { it.flatten()}

    //
    // Check mate
    //
    SMASH (
        ch_smash
    )

    emit:
    best_guesses      = SMASH.out.best_guesses     // path
}
