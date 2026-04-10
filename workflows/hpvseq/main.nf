/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SRATOOLS_FASTERQDUMP   } from '../../modules/nf-core/sratools/fasterqdump'
include { REFORMAT_FASTQ         } from '../../modules/local/reformatfastq'
include { FASTQC                 } from '../../modules/nf-core/fastqc'
include { MULTIQC                } from '../../modules/nf-core/multiqc'
include { TAGTOHEADER            } from '../../modules/local/tagtoheader'
include { GATK4_BASERECALIBRATOR } from '../../modules/nf-core/gatk4/baserecalibrator'
include { GATK4_APPLYBQSR        } from '../../modules/nf-core/gatk4/applybqsr'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_hpvseq_pipeline'

include { ALIGN_BWA              } from '../../subworkflows/local/align_bwa'
include { GATK4_BQSR             } from '../../subworkflows/local/gatk4_bqsr'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HPVSEQ {

    take:
    ch_samplesheet       // channel: samplesheet read in from --input
    ch_bwa_index         // channel: path(bwa_index/) for alignment 
    genome               // string: reference genome name, e.g. hg19 
    ch_fasta             // channel: path(genome.fasta)
    ch_fai               // channel: path(genome.fai)
    ch_dict              // channel: path(genome.dict)
    ch_bed               // channel: path(bed)
    ch_known_sites       // channel: path(known_sites)
    ch_known_sites_tbi   // channel: path(known_sites_tbi)

    main:

    ch_versions = channel.empty()
   // ch_multiqc_files = channel.empty()

    ch_samplesheet
        .branch {
            sra:   it[0].is_sra
            fastq: !it[0].is_sra
        }
        .set { ch_inputs }
    // Use Channel.value() instead of just file()
    ch_blist = params.blist ? Channel.value(file(params.blist)) : Channel.empty()


    //
    // MODULE: Run sratools/fasterdump
    //
    // Define paths (or leave empty if not needed)
    ch_ncbi_settings = [] 
    ch_certificate   = []
    SRATOOLS_FASTERQDUMP (
        ch_inputs.sra,    // tuple val(meta), path(sra)
        ch_ncbi_settings,  // path ncbi_settings
        ch_certificate     // path certificate
    )

    //
    // MODULE: reformat fastq
    //
    REFORMAT_FASTQ (
        SRATOOLS_FASTERQDUMP.out.reads
    )
    ch_all_fastqs = REFORMAT_FASTQ.out.reads.mix(ch_inputs.fastq).map { meta, r1, r2 -> [ meta, [r1, r2] ] }
    
    //
    // MODULE: Run tag_to_header  
    //
    TAGTOHEADER (
        ch_all_fastqs,
        ch_blist
    )
    ch_tag2header_reads = TAGTOHEADER.out 
    .map { meta, r1, r2 -> [ meta, [r1, r2 ] ] }
 
    //
    // MODULE: Run bwa mem  
    //
/*
    ch_fasta_fai = ch_fasta.map { fasta ->
        def fai = "${fasta}.fai" // Concatenates and turns it into a File object
        return [ [id: genome], file(fasta), file(fai) ]
    }
*/
    ch_fasta_fai          = ch_fasta.join(ch_fai)
    ALIGN_BWA (
        ch_tag2header_reads,
        ch_bwa_index, 
        ch_fasta_fai
    )

    //
    // MODULE: Run gatk4 baserecalibrator  
    //

    ch_bqsr = ALIGN_BWA.out.bam
        .join(ALIGN_BWA.out.bai)
        .combine(ch_bed)
    GATK4_BQSR (
        ch_bqsr,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_known_sites,
        ch_known_sites_tbi
    )
    emit:
    bam            = GATK4_BQSR.out.bam
    bai            = GATK4_BQSR.out.bai
    multiqc_report = "Later"
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
