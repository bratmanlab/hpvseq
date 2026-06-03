//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools'

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    input // channel: [ val(meta), [ bam ], path(fasta), path(fai) ]
    //ch_bam // channel: [ val(meta), [ bam ], path(fasta), path(fai) ]
    //ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]

    main:
    
    //SAMTOOLS_SORT(ch_bam, ch_fasta_fai, '')
    SAMTOOLS_SORT(input, '')

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.index, by: [0])
        .set { ch_bam_bai }

    ch_fasta_fai = input.map { meta, bam, fasta, fai -> tuple( meta, fasta, fai ) }    
    BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta_fai)

    emit:
    bam      = SAMTOOLS_SORT.out.bam // channel: [ val(meta), [ bam ] ]
    index    = SAMTOOLS_INDEX.out.index // channel: [ val(meta), [ index ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}
