#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/hpvseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/hpvseq
    Website: https://nf-co.re/hpvseq
    Slack  : https://nfcore.slack.com/channels/hpvseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HPVSEQ  }                 from './workflows/hpvseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hpvseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hpvseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_hpvseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta        = getGenomeAttribute('fasta')
params.fai          = getGenomeAttribute('fai')
params.dict         = getGenomeAttribute('dict')
params.bwa_index    = getGenomeAttribute('bwa')
params.dbsnp        = getGenomeAttribute('dbsnp')
params.dbsnp_tbi    = getGenomeAttribute('dbsnp_tbi')
params.indels       = getGenomeAttribute('indels')
params.indels_tbi   = getGenomeAttribute('indels_tbi')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_HPVSEQ {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HPVSEQ (
        samplesheet
    )
    // emit:
    // multiqc_report = HPVSEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
}
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.skip_fastqc,
        params.skip_multiqc,
        params.blist,
        params.genome,
        params.bed,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //

    file_fasta            = file(params.fasta, checkIfExists: true)
    file_fai              = file(params.fai, checkIfExists: true)
    file_dict             = file(params.dict, checkIfExists: true)
    ch_fasta              = channel.value([ [id: params.genome], file_fasta ])
    ch_fai                = channel.value([ [id: params.genome], file_fai ])
    ch_dict               = channel.value([ [id: params.genome], file_dict ])
    ch_bwa_index          = channel.value(file(params.bwa_index, checkIfExists: true))
    ch_bed                = channel.value(file(params.bed, checkIfExists: true))
/*
    // USE THIS INSTEAD TO LOCATE THE NULL:
    def check_params = [
        "fasta": params.fasta,
        "fai": params.fai,
        "dct": params.dict,
        "dbsnp": params.dbsnp,
        "indels": params.indels,
        "genome": params.genome,
        "bed": params.bed
    ]

    check_params.each { name, value ->
        if (value == null) {
            println "CRITICAL ERROR: parameter 'params.${name}' is NULL"
        } else {
            println "OK: params.${name} is [${value}]"
        }
    }

    println "DEBUG: dbsnp is [${params.dbsnp}]"
    println "DEBUG: indels is [${params.indels}]"
    println "DEBUG: genome is [${params.genome}]"
*/

    def vcf_paths = [params.dbsnp] + (params.indels ? params.indels.split(',') : []).flatten()
    def all_vcf = vcf_paths
        .collect { path -> file(path, checkIfExists: true) }
    ch_known_sites = channel.value([ [id: params.genome], all_vcf ])

    def tbi_paths = [params.dbsnp_tbi] + (params.indels_tbi ? params.indels_tbi.split(',') : []).flatten()
    def all_tbi = tbi_paths
        .collect { path -> file(path, checkIfExists: true) }
    ch_known_sites_tbi = channel.value([ [id: params.genome], all_tbi ])

    HPVSEQ (
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_bwa_index,
        params.genome,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_bed,
        ch_known_sites,
        ch_known_sites_tbi
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        HPVSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
