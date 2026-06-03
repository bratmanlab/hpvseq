//
// Coverage for quantification on properly-paired reads
//
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PROPERLYPAIRED     } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index'
include { GATK3_DEPTHOFCOVERAGE as COVERAGE_QUANTIFICATION_CORRECTED  } from '../../../modules/local/gatk3/depthofcoverage'

workflow COVERAGE_QUANTIFICATION_CORRECTED_F2 {
    take:
    input           // channel: [ val(meta), path(bam), path(bai), path(fasta), path(fai), path(dict), path(bed)
    //fasta_fai       // channel: [ val(meta), path(fasta), path(bai) ]
    //dict            // channel: [ val(meta), path(dict) ]
    //bed            // channel: [ val(meta), path(bed) ]

    main:

    //
    // Extract properly-paired reads (-f 2)
    //
    ch_bam = input.map { meta, bam, bai, fasta, fai, dict, bed -> tuple( meta, bam, bai ) }
    SAMTOOLS_VIEW_PROPERLYPAIRED (
        ch_bam,
        [[:], [], []],
        [[:], []],
        [[:], []],
        ""
    )
    ch_quantification_properlypaired = SAMTOOLS_VIEW_PROPERLYPAIRED.out.bam
    SAMTOOLS_INDEX(ch_quantification_properlypaired)
    ch_quantification_coverage = ch_quantification_properlypaired
    .join(SAMTOOLS_INDEX.out.index)

    ch_info = input.map { meta, bam, bai, fasta, fai, dict, bed -> tuple( meta, fasta, fai, dict, bed ) }

    ch_quantification_coverage = ch_quantification_coverage
   .join(ch_info)


    COVERAGE_QUANTIFICATION_CORRECTED (
        ch_quantification_coverage
//        ch_fasta_fai,
//        ch_dict,
//        ch_bed
    )

    emit:
    summary      = COVERAGE_QUANTIFICATION_CORRECTED.out.summary 
}

