//
// Coverage for quantification on properly-paired reads
//
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PROPERLYPAIRED     } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index'
include { GATK3_DEPTHOFCOVERAGE as COVERAGE_QUANTIFICATION_DOMINANT  } from '../../../modules/local/gatk3/depthofcoverage'

workflow COVERAGE_QUANTIFICATION_DOMINANT_F2 {
    take:
    input           // channel: [ val(meta), path(bam), path(bai) ]
    fasta_fai       // channel: [ val(meta), path(fasta), path(bai) ]
    dict            // channel: [ val(meta), path(dict) ]
    bed            // channel: [ val(meta), path(bed) ]

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
    ch_quantification_properlypaired = SAMTOOLS_VIEW_PROPERLYPAIRED.out.bam
    SAMTOOLS_INDEX(ch_quantification_properlypaired)
    ch_quantification_coverage = ch_quantification_properlypaired.join(SAMTOOLS_INDEX.out.index)

    COVERAGE_QUANTIFICATION_DOMINANT (
        ch_quantification_coverage,
        fasta_fai,
        dict,
        bed
    )

    emit:
    summary      = COVERAGE_QUANTIFICATION_DOMINANT.out.summary 
}

