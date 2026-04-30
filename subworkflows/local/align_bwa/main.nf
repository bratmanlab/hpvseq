//
// Alignment with BWA MEM
//

include { BWA_MEM           } from '../../../modules/nf-core/bwa/mem'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools'

//
// Function that parses and returns the alignment rate from the BWA MEM log output
//
def getBwaPercentMapped(align_flagstat) {
    def percent_mapped = 0
    // def line = "4534962 + 0 mapped (61.14% : N/A)"
    def pattern = /mapped \(([\d.]+)%/
    align_flagstat.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_mapped = matcher[0][1].toFloat()
        }
    }
    return percent_mapped
}

workflow ALIGN_BWA {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    index         // channel: /path/to/bwa/index/ ?
    fasta         // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    //
    // Map reads with BWA MEM
    //
    BWA_MEM(
        reads,
        index.map { index_path -> [ [id: 'genome'], index_path ] },
        [ [:], [] ],    // No fasta needed for BAM output
        false           // sort_bam - we'll sort with samtools for consistency
    )
    ch_orig_bam = BWA_MEM.out.bam

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    //BAM_SORT_STATS_SAMTOOLS(ch_orig_bam, fasta)
    BAM_SORT_STATS_SAMTOOLS(
        ch_orig_bam, 
        ch_orig_bam
        .combine(fasta)
        .map { meta_id, bam, meta_genome, fasta, fai ->
            [ meta_id, fasta, fai]
        }.first() 
    ) 
    ch_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat

    //
    // Parse alignment rate from log
    //
    ch_percent_mapped = ch_flagstat.map { meta, log_file -> [ meta, getBwaPercentMapped(log_file) ] }

    emit:
//    orig_bam       = ch_orig_bam                          // channel: [ val(meta), bam ]
    bam            = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai            = BAM_SORT_STATS_SAMTOOLS.out.index      // channel: [ val(meta), [ bai ] ]
    stats          = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat       = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats       = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    percent_mapped = ch_percent_mapped                    // channel: [ val(meta), percent_mapped ]
}
