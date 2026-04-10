//
// GATK4 BQSR: BASERECALIBRATOR and APPLYBQSR
//

include { GATK4_BASERECALIBRATOR  } from '../../../modules/nf-core/gatk4/baserecalibrator'
include { GATK4_APPLYBQSR         } from '../../../modules/nf-core/gatk4/applybqsr'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools'

workflow GATK4_BQSR {
    take:
    bam_interval         // channel: [ val(meta), path(input), path(input_index), path(intervals) ]
    fasta         // channel: /path/to/fasta
    fai         // channel: /path/to/fai
    dict         // channel: /path/to/dict
    known_sites         // channel: /path/to/known_sites
    known_sites_tbi         // channel: /path/to/known_sites_tbi

    main:

    //
    // Baserecalibrator
    //
    GATK4_BASERECALIBRATOR (
        bam_interval,
        fasta,
        fai,
        dict,
        known_sites,
        known_sites_tbi
    )
    
    //
    // ApplyBQSR without interval to keep full reads 
    //
    ch_bam = bam_interval
        .join(GATK4_BASERECALIBRATOR.out.table)
        .map { meta, bam, bai, intervals, table ->
            [ meta, bam, bai, table, [] ]
    }

    file_fasta = fasta.map { meta, path -> path }
    file_fai = fai.map { meta, path -> path }
    file_dict = dict.map { meta, path -> path }
    GATK4_APPLYBQSR (
        ch_bam,
        file_fasta,
        file_fai,
        file_dict
    )

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_fasta_fai = fasta.join(fai)
    BAM_SORT_STATS_SAMTOOLS(GATK4_APPLYBQSR.out.bam, ch_fasta_fai)
    ch_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat

    emit:
    bam            = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai            = BAM_SORT_STATS_SAMTOOLS.out.index      // channel: [ val(meta), [ bai ] ]
    stats          = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat       = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
}
