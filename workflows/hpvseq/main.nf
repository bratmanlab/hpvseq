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
include { CONSENSUSCRUNCHER      } from '../../modules/local/consensuscruncher'
include { PICARD_COLLECTHSMETRICS} from '../../modules/nf-core/picard/collecthsmetrics'
include { GATK3_DEPTHOFCOVERAGE as COVERAGE_QC_BWA          } from '../../modules/local/gatk3/depthofcoverage'
include { GATK3_DEPTHOFCOVERAGE as COVERAGE_QC_DCS          } from '../../modules/local/gatk3/depthofcoverage'
include { GATK3_DEPTHOFCOVERAGE as COVERAGE_QC_DCSSC        } from '../../modules/local/gatk3/depthofcoverage'
include { GATK3_DEPTHOFCOVERAGE as COVERAGE_QC_SSCS         } from '../../modules/local/gatk3/depthofcoverage'
include { GATK3_DEPTHOFCOVERAGE as COVERAGE_QC_SSCSSC       } from '../../modules/local/gatk3/depthofcoverage'
include { GATK3_DEPTHOFCOVERAGE as COVERAGE_QC_ALLUNIQUEDCS } from '../../modules/local/gatk3/depthofcoverage'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PROPERLYPAIRED     } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX         } from '../../modules/nf-core/samtools/index'
include { SUMMARIZE_REPORT       } from '../../modules/local/summarizereport'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_hpvseq_pipeline'

include { ALIGN_BWA              } from '../../subworkflows/local/align_bwa'
include { GATK4_BQSR             } from '../../subworkflows/local/gatk4_bqsr'
include { GENOTYPING             } from '../../subworkflows/local/genotyping'
include { QUANTIFICATION         } from '../../subworkflows/local/quantification'
include { COVERAGE_QUANTIFICATION_F2 as COVERAGE_QUANTIFICATION_HG    } from '../../subworkflows/local/coverage_quantification_f2'
include { COVERAGE_QUANTIFICATION_F2 as COVERAGE_QUANTIFICATION_VIRUS } from '../../subworkflows/local/coverage_quantification_f2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HPVSEQ {

    take:
    ch_samplesheet       // channel: samplesheet read in from --input
    ch_bwa_index         // channel: path(bwa_index/) for alignment 
    genome               //  string: reference genome name, e.g. hg19 
    ch_fasta             // channel: path(genome.fasta)
    ch_fai               // channel: path(genome.fai)
    ch_dict              // channel: path(genome.dict)
    ch_bed               // channel: path(bed)
    sequencing_platform  //  string: sequencing platform
    consensuscruncher_dir  //  string: path(consensuscruncher) 
    ch_known_sites       // channel: path(known_sites)
    ch_known_sites_tbi   // channel: path(known_sites_tbi)
    ch_genotyping_index  // channel: path(genotyping index) for alignment for genotyping
    ch_genotyping_fasta  // channel: path(genotyping fasta) for alignment for genotyping
    ch_genotyping_fai    // channel: path(genotyping fai) for alignment for genotyping
    ch_genotyping_cytoband    // channel: path(genotyping cytoband) for alignment for genotyping
    ch_genotypes         // channel: path(genotypes) for alignment for detected genotype

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
    ch_blist = params.blist ? channel.value(file(params.blist)) : channel.empty()


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

    // Add RG to meta
    ch_all_fastqs_rg = ch_all_fastqs.map { meta, reads ->
        def fastq = reads[0] // Check R1
        // Extract first line: @Instrument:Run:Flowcell:Lane:Tile:X:Y Read:Filter:Control:Index
        def header = "zcat ${fastq} | head -n 1".execute().text.trim()
        
        // Parse the sequencer header (e.g. Illumina) (adjust regex if using non-standard headers)
        def parts = header.split(':')
        if (parts.size() >= 7) {
            def flowcell = parts[2]
            def lane     = parts[3]
            def index = "unknown"
            if (parts.size() >= 10) {//e.g. from SRR: M05097:166:000000000-KYHTV:1:1101:15530:1563 length=151
                index    = parts[9] // Handle space before second part of header
            }
            //rg=$(echo "@RG\tID:"$flowcell"_"$lane"\tSM:"$samp"\tPL:Illumina\tPU:.\tLB:"$lib)
            meta = meta + [ rg: "\"@RG\\tID:${flowcell}_${lane}\\tSM:${meta.id}\\tPL:${sequencing_platform}\\tPU:${flowcell}.${lane}.${index}\\tLB:${index}\"" ]
        }else{
            meta = meta + [ rg: "\"@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:${sequencing_platform}\\tPU:.\\tLB:.\"" ]
        }
        return [ meta, reads ]
    }
    
    //
    // MODULE: Run fastqc and multiqc  
    //
    FASTQC (
        ch_all_fastqs_rg
    )

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = FASTQC.out.zip
    .map { meta, file -> file } // Explicitly keep only the file, discard meta
    .collect()  
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] })
    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    // MULTIQC 
    MULTIQC (
        ch_multiqc_files, // Aggregate all files into one list
        [],                         // multiqc_config (optional)
        [],                         // extra_multiqc_config (optional)
        [],                         // multiqc_logo (optional)
        [],                         // multiqc_replace_names (optional)
        []                          // multiqc_title (optional)
    )

    //
    // MODULE: Run merge fastq files 
    //
    CAT_FASTQ (
        ch_all_fastqs_rg
    )

    //
    // MODULE: Run trim read length 
    //


    //
    // MODULE: Run tag_to_header  
    //
    TAGTOHEADER (
        ch_all_fastqs_rg,
        ch_blist
    )

    ch_fasta_fai          = ch_fasta.join(ch_fai)
    ch_tag2header_reads = TAGTOHEADER.out
    .combine(ch_fasta_fai) 
    .map { meta, r1, r2, meta_genome, fasta, fai -> 
        meta = meta + [ref: meta_genome.id, type: "Aligned", ref2: "virus", type2: "Unknown", consensus: "none" ]
        return [ meta, [r1, r2 ] ] 
    }
 
    //
    // MODULE: Run bwa mem  
    //
    ch_fasta_fai          = ch_fasta.join(ch_fai)
    ALIGN_BWA (
        ch_tag2header_reads,
        ch_bwa_index, 
        ch_fasta_fai
    )
    ch_bam = ALIGN_BWA.out.bam
    .join(ALIGN_BWA.out.bai)
/*
    ch_bam.view { it -> 
        println(it)
        println(it*.getClass())
    }
*/
    // Report: on-target rate 
    ch_hsmetrics = ch_bam
    .map { meta, bam, bai ->
        return [ meta, bam, bai, file(params.bait_intervals, checkIfExists: true), file(params.target_intervals, checkIfExists: true) ]
    }

    PICARD_COLLECTHSMETRICS (
        ch_hsmetrics,
        ch_fasta,
        ch_fai,
        ch_dict,
        [ [:], [] ]    
    ) 

    // Report: QC coverage (without -f 2) 
    COVERAGE_QC_BWA (
        ch_bam,
        ch_fasta_fai.first(),
        ch_dict,
        ch_bed
    )
 
    //
    // MODULE: Run gatk4 baserecalibrator  
    //

    ch_bqsr = ch_bam
    .map { meta, bam, bai ->
        // Use params.intervals if it exists, otherwise use an empty list []
        // This ensures the 4th element is NEVER null
        def intervals = params.bed ? file(params.bed) : []
        return [ meta, bam, bai, intervals ]
    }

    GATK4_BQSR (
        ch_bqsr,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_known_sites,
        ch_known_sites_tbi
    )
    
    //
    // ConsensusCruncher:HG
    //
    ch_cc_hg = GATK4_BQSR.out.bam
       .join(GATK4_BQSR.out.bai)
    file_cytoband = file("${params.consensuscruncher_dir}/ConsensusCruncher/${params.genome}_cytoBand.txt", checkIfExists: true)
    ch_cytoband = channel.value( [ [id: params.genome], file_cytoband ] )
/*
    ch_cc_hg.view { meta, bam, bai -> 
	"""
	Sample ID: ${meta.id}
	Reference: ${meta.ref}
	BAM type : ${meta.type}
	BAM Path : ${bam}
	BAI Path : ${bai}
	-----------------
	""".stripIndent()
    }
*/    
    CONSENSUSCRUNCHER (
        ch_cc_hg,
        ch_cytoband
    )
    ch_dcs           = CONSENSUSCRUNCHER.out.bam_dcs.join(CONSENSUSCRUNCHER.out.bai_dcs)
    ch_dcssc         = CONSENSUSCRUNCHER.out.bam_dcssc.join(CONSENSUSCRUNCHER.out.bai_dcssc)
    ch_sscs          = CONSENSUSCRUNCHER.out.bam_sscs.join(CONSENSUSCRUNCHER.out.bai_sscs)
    ch_sscssc        = CONSENSUSCRUNCHER.out.bam_sscssc.join(CONSENSUSCRUNCHER.out.bai_sscssc)
    ch_alluniquedcs  = CONSENSUSCRUNCHER.out.bam_alluniquedcs.join(CONSENSUSCRUNCHER.out.bai_alluniquedcs)
    ch_dcs           = ch_dcs.map { meta, bam, bai -> [ meta + [ consensus: "dcs" ], bam ,bai ] }
    ch_dcssc         = ch_dcssc.map { meta, bam, bai -> [ meta + [ consensus: "dcs.sc" ], bam ,bai ] }
    ch_sscs          = ch_sscs.map { meta, bam, bai -> [ meta + [ consensus: "sscs" ], bam ,bai ] }
    ch_sscssc        = ch_sscssc.map { meta, bam, bai -> [ meta + [ consensus: "sscs.sc" ], bam ,bai ] }
    ch_alluniquedcs  = ch_alluniquedcs.map { meta, bam, bai -> [ meta + [ consensus: "all.unique.dcs" ], bam ,bai ] }

    // Report: QC coverage (without -f 2)
    COVERAGE_QC_DCS ( ch_dcs, ch_fasta_fai.first(), ch_dict, ch_bed)
    COVERAGE_QC_DCSSC ( ch_dcssc, ch_fasta_fai.first(), ch_dict, ch_bed)
    COVERAGE_QC_SSCS ( ch_sscs, ch_fasta_fai.first(), ch_dict, ch_bed)
    COVERAGE_QC_SSCSSC ( ch_sscssc, ch_fasta_fai.first(), ch_dict, ch_bed)
    COVERAGE_QC_ALLUNIQUEDCS ( ch_alluniquedcs, ch_fasta_fai.first(), ch_dict, ch_bed)

    // Report: coverage quantification hg
    COVERAGE_QUANTIFICATION_HG (
        ch_alluniquedcs,
        ch_fasta_fai.first(),
        ch_dict,
        ch_bed 
    )

    //
    // MODULE: Run genotyping 
    //
    ch_genotyping_fasta_fai  = ch_genotyping_fasta.join(ch_genotyping_fai)

    GENOTYPING (
        ch_bam,
        ch_genotyping_index, 
        ch_genotyping_fasta_fai,
        ch_genotyping_cytoband,
        ch_genotypes
    )

    ch_quantification = GENOTYPING.out.quantification_input

    //
    // QUANTIFICATION: re-align unmapped reads on given genotype
    //
    QUANTIFICATION (
        ch_quantification
    )

    ch_quantification_reads = QUANTIFICATION.out.bam_alluniquedcs.join(QUANTIFICATION.out.bai_alluniquedcs)
    ch_quantification_reads = ch_quantification_reads
    .map { meta, bam, bai -> [ meta + [ consensus: "all.unique.dcs"], bam, bai ] }

    // Report: Quantification coverage (with -f 2)
    COVERAGE_QUANTIFICATION_VIRUS (
        ch_quantification_reads,
        QUANTIFICATION.out.fasta_fai_genotype,
        QUANTIFICATION.out.dict_genotype,
        QUANTIFICATION.out.bed_genotype
    )
 
    //
    // Summary: on-target rate (hsmetrics), coverage qc, genotyping, coverage quantification, saturation rate
    //

    ch_hsmetrics_report = PICARD_COLLECTHSMETRICS.out.metrics.map { meta, file -> file }.collect()

    ch_coverage_qc_bwa_report           = COVERAGE_QC_BWA.out.summary.collect()
/*
    ch_coverage_qc_bwa_report.view { it ->
        println(it)
        println(it*.getClass())
    }
*/
    ch_coverage_qc_dcs_reprot           = COVERAGE_QC_DCS.out.summary.collect()
    ch_coverage_qc_dcssc_report         = COVERAGE_QC_DCSSC.out.summary.collect()
    ch_coverage_qc_sscs_report          = COVERAGE_QC_SSCS.out.summary.collect()
    ch_coverage_qc_sscssc_report        = COVERAGE_QC_SSCSSC.out.summary.collect()
    ch_coverage_qc_alluniquedcs_report  = COVERAGE_QC_ALLUNIQUEDCS.out.summary.collect()

    ch_genotyping_report = GENOTYPING.out.genotyping_reads.map { meta, file -> file }.collect()
    genotypes_info       = ch_genotypes.map { map, file -> file}

    ch_coverage_quantification_hg       = COVERAGE_QUANTIFICATION_HG.out.summary.collect()
    ch_coverage_quantification_virus    = COVERAGE_QUANTIFICATION_VIRUS.out.summary.collect()

    ch_read_families        = CONSENSUSCRUNCHER.out.optional_read_families.collect()
    ch_read_families_virus  = QUANTIFICATION.out.read_families.collect()

    SUMMARIZE_REPORT (
        ch_hsmetrics_report,
        ch_coverage_qc_bwa_report,
        ch_coverage_qc_dcs_reprot,
        ch_coverage_qc_dcssc_report,
	ch_coverage_qc_sscs_report,
	ch_coverage_qc_sscssc_report,
	ch_coverage_qc_alluniquedcs_report,
	ch_genotyping_report,
        genotypes_info,
        ch_coverage_quantification_hg, 
        ch_coverage_quantification_virus,
        ch_read_families, 
        ch_read_families_virus 
    )     

    emit:
    cc_dir                  = CONSENSUSCRUNCHER.out.consensus_dir
    cc_dir_genotype         = QUANTIFICATION.out.consensus_dir
    multiqc_report          = "Later"
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
