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

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_hpvseq_pipeline'

include { ALIGN_BWA              } from '../../subworkflows/local/align_bwa'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HPVSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_bwa_index      // channel: path(bwa_index/) for alignment 
    ch_fasta       // channel: path(genome.fasta)

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

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
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
 
    //
    // MODULE: Run tag_to_header  
    //
    TAGTOHEADER (
        ch_all_fastqs,
        ch_blist
    )
   
    //
    // MODULE: Run bwa mem  
    //
    ALIGN_BWA (
        TAGTOHEADER.out.reads,
        ch_bwa_index, 
        ch_fasta.map { item -> [ [:], item ] }
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
