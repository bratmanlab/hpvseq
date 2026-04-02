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
    // def pattern = /(\d+\.\d+)% overall alignment rate/ // Bowtie2
    def pattern = /(\d+\.\d+)% mapped/
    align_flagstat.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_mapped = matcher[0][1].toFloat()
        }
    }
    return percent_aligned
}

workflow ALIGN_BWA {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    index         // channel: /path/to/bwa/index/
    fasta         // channel: /path/to/fasta

    main:

    //
    // Map reads with BWA MEM
    //
    BWA_ALIGN(
        reads,
        index.map { index_path -> [ [id: 'genome'], index_path ] },
        [ [:], [] ],    // No fasta needed for BAM output
        false           // sort_bam - we'll sort with samtools for consistency
    )

    ch_orig_bam = BWA_MEM.out.bam
    ch_log = BWA_MEM.out.log

    // Parse alignment rate from log
    ch_percent_mapped = ch_log.map { meta, log_file -> [ meta, getBwaPercentMapped(log_file) ] }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(ch_orig_bam, fasta)
    ch_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat

    //
    // Parse alignment rate from log
    //
    ch_percent_mapped = ch_flagstat.map { meta, log_file -> [ meta, getBwaPercentMapped(log_file) ] }

    emit:
    orig_bam       = ch_orig_bam                          // channel: [ val(meta), bam ]
    log_final      = ch_log                               // channel: [ val(meta), log ]
    bam            = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai            = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi            = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats          = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat       = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats       = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    percent_mapped = ch_percent_mapped                    // channel: [ val(meta), percent_mapped ]
}
