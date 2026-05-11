/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SRATOOLS_FASTERQDUMP   } from '../../modules/nf-core/sratools/fasterqdump'
include { INDEX_MISMATCH         } from '../../modules/local/index_mismatch'
include { META_TIMEPOINT         } from '../../modules/local/meta_timepoint'
include { CAT_FASTQ              } from '../../modules/nf-core/cat/fastq'
include { FASTP                  } from '../../modules/nf-core/fastp'
include { REFORMAT_FASTQ         } from '../../modules/local/reformatfastq'
include { REFORMAT_FASTQ2        } from '../../modules/local/reformatfastq2'
include { FASTQC                 } from '../../modules/nf-core/fastqc'
include { MULTIQC                } from '../../modules/nf-core/multiqc'
include { TAGTOHEADER            } from '../../modules/local/tagtoheader'
include { COUNT_READS as COUNT_READS_RAW      } from '../../modules/local/count_reads'
include { COUNT_READS as COUNT_READS_GOODUMI  } from '../../modules/local/count_reads'
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
include { INTEGRATION_SEARCHPV   } from '../../modules/local/integration/searchpv'
include { SUMMARIZE_TOHSMETRICS  } from '../../modules/local/summarize/tohsmetrics'
include { SUMMARIZE_REPORTS      } from '../../modules/local/summarize/reports'
include { SUMMARIZE_QUANTIFICATION_CORRECTED        } from '../../modules/local/summarize/quantification_corrected'
include { SUMMARIZE_INSERTSIZE   } from '../../modules/local/summarize/insertsize'
include { SUMMARIZE_INTEGRATION  } from '../../modules/local/summarize/integration'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_hpvseq_pipeline'

include { CHECKMATE              } from '../../subworkflows/local/checkmate'
include { ALIGN_BWA              } from '../../subworkflows/local/align_bwa'
include { GATK4_BQSR             } from '../../subworkflows/local/gatk4_bqsr'
include { GENOTYPING             } from '../../subworkflows/local/genotyping'
include { INPUT_QUANTIFICATION_CORRECTED             } from '../../subworkflows/local/input/quantification_corrected'
include { QUANTIFICATION as QUANTIFICATION_DOMINANT  } from '../../subworkflows/local/quantification'
include { QUANTIFICATION as QUANTIFICATION_CORRECTED } from '../../subworkflows/local/quantification'
include { COVERAGE_QUANTIFICATION_DOMINANT_F2 as COVERAGE_QUANTIFICATION_HG    } from '../../subworkflows/local/coverage_quantification_dominant_f2'
include { COVERAGE_QUANTIFICATION_DOMINANT_F2 as COVERAGE_QUANTIFICATION_VIRUS } from '../../subworkflows/local/coverage_quantification_dominant_f2'
include { COVERAGE_QUANTIFICATION_CORRECTED_F2 as COVERAGE_QUANTIFICATION_CORRECTED_VIRUS } from '../../subworkflows/local/coverage_quantification_corrected_f2'
include { INSERTSIZE as INSERTSIZE_HG         } from '../../subworkflows/local/insertsize'
include { INSERTSIZE as INSERTSIZE_DOMINANT   } from '../../subworkflows/local/insertsize'
include { INSERTSIZE as INSERTSIZE_CORRECTED  } from '../../subworkflows/local/insertsize'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HPVSEQ {

    take:
    ch_samplesheet         // channel: samplesheet read in from --input
    ch_index               // channel: path(i7_i5 index)
    ch_bwa_index           // channel: [ meta, path(bwa_index/) ] 
    genome                 //  string: reference genome name, e.g. hg19 
    ch_fasta               // channel: [ meta, path(genome.fasta) ]
    ch_fai                 // channel: [ meta, path(genome.fai) ]
    ch_dict                // channel: [ meta, path(genome.dict) ]
    ch_bed                 // channel: [ meta, path(bed) ]
    sequencing_platform    //  string: sequencing platform
    consensuscruncher_dir  //  string: path(consensuscruncher) 
    ch_known_sites         // channel: [ meta, path(known_sites) ]
    ch_known_sites_tbi     // channel: [ meta, path(known_sites_tbi) ]
    ch_genotyping_index    // channel: [ meta2, path(genotyping index) ] - for alignment for genotyping 
    ch_genotyping_fasta    // channel: [ meta2, path(genotyping fasta) ] - for alignment for genotyping 
    ch_genotyping_fai      // channel: [ meta2, path(genotyping fai) ] - for alignment for genotyping 
    ch_genotyping_cytoband    // channel: [ meta2,path(genotyping cytoband) ] - for alignment for genotyping 
    ch_genotypes           // channel: [ meta3, path(genotypes) ] - for alignment for detected genotype

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
    ch_all_fastqs = REFORMAT_FASTQ.out.reads.mix(ch_inputs.fastq).map { meta, r1, r2 -> tuple( meta, [r1, r2] ) }

    // Add RG to meta
    ch_all_fastqs_rg = ch_all_fastqs.map { meta, reads ->
        def fastq = reads[0] // Check R1
        // Extract first line: @Instrument:Run:Flowcell:Lane:Tile:X:Y Read:Filter:Control:Index
        //def header = "zcat ${fastq} | head -n 1".execute().text.trim()
        def header = ["bash", "-c", "zcat ${fastq} | head -n 1"].execute().text.trim()
        
        // Parse the sequencer header (e.g. Illumina) (adjust regex if using non-standard headers)
        def parts = header.split(':')
        if (parts.size() >= 7) {
            def flowcell = parts[2]
            //def lane     = parts[3]
            def index = "unknown"
            if (parts.size() >= 10) {//e.g. from SRR: M05097:166:000000000-KYHTV:1:1101:15530:1563 length=151
                index    = parts[9] // Handle space before second part of header
            }
            //rg=$(echo "@RG\tID:"$flowcell"_"$lane"\tSM:"$samp"\tPL:Illumina\tPU:.\tLB:"$lib)
            meta = meta + [ rg: "\"@RG\\tID:${flowcell}\\tSM:${meta.id}\\tPL:${sequencing_platform}\\tPU:${flowcell}.${index}\\tLB:${index}\"" ]
        }else{
            meta = meta + [ rg: "\"@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:${sequencing_platform}\\tPU:.\\tLB:.\"" ]
        }
        tuple( meta, reads )
    }
    
    //
    // MODULE: Check sample index mismatch  
    //
    ch_index_mismatch = ch_all_fastqs_rg
    .map { meta, reads ->
        tuple( meta + [ id: "${meta.id}_L${meta.lane}"], reads ) 
    }
    INDEX_MISMATCH (
        ch_index_mismatch,
        ch_index,
        params.index_mismatch
    )

    INDEX_MISMATCH.out.reads
        .branch { meta, reads, passfile, indexfile ->
            passed: passfile.text.trim() == "1"
            failed: passfile.text.trim() == "0"
        }
        .set { ch_status }

    // Print a warning for every failed sample
    ch_status.failed.subscribe { meta, reads, passfile, indexfile ->
        log.warn "!!! [WARNING] Sample ${meta.id} FAILED validation of INDEX MISMATCH and will be skipped."
    }

    //
    // MODULE: Update meta tp
    //
    ch_passed = ch_status.passed
    .map { meta, reads, passfile, indexfile ->
        tuple( meta, reads, indexfile ) 
    }
    META_TIMEPOINT (
        ch_passed,
        params.library
    )
    ch_fastqc = META_TIMEPOINT.out.reads
    .map { meta, reads, tpfile ->
        def tp = file(tpfile).text.trim()
        tuple( meta + [ tp: tp ], reads ) 
    }

    //
    // MODULE: Run fastqc and multiqc  
    //
    FASTQC (
        ch_fastqc
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
    ch_cat_fastq = ch_fastqc
    .map { meta, reads ->
	def new_meta = meta.clone()
        new_meta.id = meta.id.replaceAll(/_L\d+$/, "")
        new_meta.remove('lane') // VERY IMPORTANT: remove the lane key
        tuple( new_meta, reads ) 
    }
    .groupTuple()
    .map { meta, reads -> 
        // This turns [[R1, R2], [R3, R4]] into [R1, R2, R3, R4]
        tuple( meta, reads.flatten() ) 
    }
    CAT_FASTQ (
        ch_cat_fastq
    )
    ch_fastqs = CAT_FASTQ.out.reads
    .map{ meta, reads ->
        def size_gb = reads.collect { it.size() }.sum() / 1024**3
        tuple( meta + [ size_gb: size_gb ], reads )
    }
    //
    // Run or skip trim read length
    // 
    if ( !params.skip_trim_read ) {
        FASTQP (
	    ch_fastqs.map { meta, reads -> tuple( meta, reads, [] ) }, // meta, reads, adapter
	    false, // discard_trimmed_pass (keep the data)
	    false, // save_trimmed_fail (don't save bad reads)
	    false  // save_merged (don't stitch R1/R2)
        )
        ch_tag2header = FASTP.out.reads
    } else {
        ch_tag2header = ch_fastqs 
    }


    //
    // MODULE: Run tag_to_header  
    //
    TAGTOHEADER (
        ch_tag2header,
        ch_blist
    )

    ch_fasta_fai          = ch_fasta.join(ch_fai)
    ch_tag2header_reads   = TAGTOHEADER.out
    .combine(ch_fasta_fai) 
    .map { meta, r1, r2, meta_genome, fasta, fai -> 
        meta = meta + [ref: meta_genome.id, type: "Aligned", ref2: "virus", type2: "Unknown", consensus: "none" ]
        tuple( meta, [r1, r2 ] ) 
    }

    // Report #read pairs
    COUNT_READS_RAW (
        ch_tag2header
    ) 
    COUNT_READS_GOODUMI (
        ch_tag2header_reads
    ) 
    //
    // MODULE: Run bwa mem  
    //
    ALIGN_BWA (
        ch_tag2header_reads,
        ch_bwa_index, 
        ch_fasta_fai
    )
    ch_bam = ALIGN_BWA.out.bam
    .join(ALIGN_BWA.out.bai)

    // Report: on-target rate 
    ch_hsmetrics = ch_bam
    .map { meta, bam, bai ->
        tuple( meta, bam, bai, file(params.bait_intervals, checkIfExists: true), file(params.target_intervals, checkIfExists: true) )
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

    // Summary: if skip_genotyping
    ch_nreads_raw_report     = COUNT_READS_RAW.out.nreads.map { meta, file -> file }.collect()
    ch_nreads_goodumi_report = COUNT_READS_GOODUMI.out.nreads.map { meta, file -> file }.collect()
    ch_hsmetrics_report      = PICARD_COLLECTHSMETRICS.out.metrics.map { meta, file -> file }.collect()
    println ch_nreads_raw_report 
    println ch_nreads_goodumi_report
 
    SUMMARIZE_TOHSMETRICS (
	ch_nreads_raw_report,
	ch_nreads_goodumi_report,
	ch_hsmetrics_report
    )     

    //
    // MODULE: Run sample swap based on de-dup bam files  
    //
    ch_checkmate_best_guesses = channel.empty()
    if (!params.skip_checkmate) {
	ch_checkmate = ch_bam.map { meta, bam, bai -> tuple( meta, bam ) }
	CHECKMATE (
	    ch_checkmate,
	    ch_fasta_fai.first()
	)
        ch_checkmate_best_guesses = CHECKMATE.out.best_guesses
    } else if (!params.skip_genotype_correction) {
        ch_checkmate_best_guesses = channel.value( file( params.checkmate_best_guesses, checkIfExists: true ) )
    }
 
    //
    // MODULE: Run gatk4 baserecalibrator  
    //

    ch_bqsr = ch_bam
    .map { meta, bam, bai ->
        // Use params.intervals if it exists, otherwise use an empty list []
        // This ensures the 4th element is NEVER null
        def intervals = params.bed ? file(params.bed) : []
        tuple( meta, bam, bai, intervals )
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
    ch_dcs           = ch_dcs.map { meta, bam, bai -> tuple( meta + [ consensus: "dcs" ], bam ,bai ) }
    ch_dcssc         = ch_dcssc.map { meta, bam, bai -> tuple( meta + [ consensus: "dcs.sc" ], bam ,bai ) }
    ch_sscs          = ch_sscs.map { meta, bam, bai -> tuple( meta + [ consensus: "sscs" ], bam ,bai ) }
    ch_sscssc        = ch_sscssc.map { meta, bam, bai -> tuple( meta + [ consensus: "sscs.sc" ], bam ,bai ) }
    ch_alluniquedcs  = ch_alluniquedcs.map { meta, bam, bai -> tuple( meta + [ consensus: "all.unique.dcs" ], bam ,bai ) }

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

    // Report: fragment insert size hg
    ch_insertsize_report_hg   = channel.empty()
    if ( !params.skip_insertsize ){
        ch_insertsize_hg = ch_alluniquedcs
        .map { meta, bam, bai ->
            tuple( meta + [ insertsize: "hg" ], bam, bai )
        }
        INSERTSIZE_HG (
            ch_insertsize_hg
        )
        ch_insertsize_report_hg  = INSERTSIZE_HG.out.metrics.map { meta, file -> file }.collect()
    }

    //
    // MODULE: Run genotyping 
    //
    ch_genotyping_fasta_fai              = ch_genotyping_fasta.join(ch_genotyping_fai)
    ch_quantification                    = channel.empty()
    ch_quantification_corrected_unmapped = channel.empty() 
    ch_genotyping_report                 = channel.empty()
    if (!params.skip_genotyping) {
	GENOTYPING (
	    ch_bam,
	    ch_genotyping_index, 
	    ch_genotyping_fasta_fai,
	    ch_genotyping_cytoband,
	    ch_genotypes
	)

	ch_quantification           = GENOTYPING.out.quantification_input
        ch_quantification_corrected_unmapped = GENOTYPING.out.bam_unmapped
        ch_genotyping_report        = GENOTYPING.out.genotyping_reads.map { meta, file -> file }.collect()
    }

    //
    // MODULE: Run SearcHPV for HPV integration on diminant genotypes
    //
    ch_integration_report            = channel.empty()
    if (!params.skip_integration && params.skip_genotype_correction) {
        REFORMAT_FASTQ2 (
            ch_tag2header_reads
        )
        ch_tag2header_reads_v2 = REFORMAT_FASTQ2.out.reads 
        .map { meta, r1, r2 ->
            def cleaned_meta = meta.subMap(['id'])
            tuple( cleaned_meta, [ r1, r2 ] ) 
        }
        ch_quantification_v2 = ch_quantification
        .map { meta, bam_unmapped, genotype_file ->
            def cleaned_meta = meta.subMap(['id'])
            tuple( cleaned_meta, genotype_file )
        }
        ch_integration_input = ch_quantification_v2
        .join(ch_tag2header_reads_v2)
        .map { meta, genotype_file, reads ->
            tuple( meta, reads, genotype_file ) 
        }
	ch_integration_input.view { it ->
            println(it)
            println(it*.getClass())
	}
        
        ch_genotype_bwa_index = ch_integration_input
        .map { meta, reads, genotypef ->
            def virus_name = ""
            if (genotypef instanceof Path || genotypef instanceof File) {
                def lines = genotypef.readLines()
                virus_name = lines ? lines[0].trim() : null
            } else {
                virus_name = genotypef.toString().trim()
            }
            def index = "${params.ref_path_virus}/${virus_name}/"
            tuple( meta + [ virus: virus_name ], file(index, checkIfExists: true) ) 
        }
        ch_integration_input = ch_integration_input.map { meta, reads, genotypef -> tuple( meta, reads ) }
        INTEGRATION_SEARCHPV (
            ch_integration_input,
            ch_bwa_index,
            ch_genotype_bwa_index 
        )
        ch_integration_report = INTEGRATION_SEARCHPV.out.HPVfusionPointContigSrNum.collect() 
    }

    //
    // QUANTIFICATION: re-align unmapped reads on given genotype
    //
    QUANTIFICATION_DOMINANT (
        ch_quantification
    )

    ch_quantification_reads = QUANTIFICATION_DOMINANT.out.bam_alluniquedcs.join(QUANTIFICATION_DOMINANT.out.bai_alluniquedcs)
    ch_quantification_reads = ch_quantification_reads
    .map { meta, bam, bai -> tuple( meta + [ consensus: "all.unique.dcs"], bam, bai ) }

    // Report: Quantification coverage (with -f 2)
    COVERAGE_QUANTIFICATION_VIRUS (
        ch_quantification_reads,
        QUANTIFICATION_DOMINANT.out.fasta_fai_genotype,
        QUANTIFICATION_DOMINANT.out.dict_genotype,
        QUANTIFICATION_DOMINANT.out.bed_genotype
    )
 
    // Report: fragment insert size virus on dominant genotype
    ch_insertsize_report_dominant_virus   = channel.empty()
    if ( !params.skip_insertsize ){
        ch_insertsize_dominant = ch_quantification_reads
        .map { meta, bam, bai ->
            tuple( meta + [ insertsize: "dominant" ], bam, bai )
        }
        INSERTSIZE_DOMINANT (
            ch_insertsize_dominant
        )
        ch_insertsize_report_dominant_virus  = INSERTSIZE_DOMINANT.out.metrics.map { meta, file -> file }.collect()
    }

    //
    // QUANTIFICATION: re-align unmapped reads on given baseline-corrected genotype 
    //
    ch_coverage_quantification_corrected_virus = channel.empty()
    ch_read_families_corrected_virus           = channel.empty()
    ch_insertsize_report_corrected_virus       = channel.empty()
    if ( !params.skip_genotype_correction ) {
        ch_quantification
        .branch { meta, bam_unmapped, genotype_file ->
            baseline : meta.tp == params.baseline_name
            others   : meta.tp != params.baseline_name
        }
        .set { ch_genotype_correction }
        ch_baseline_info = ch_genotype_correction.baseline
        .map { meta, bam_unmapped, genotype_file ->
            tuple( meta.id, genotype_file.text.trim() ) 
        }
        .toList()
        .map { list ->
            println list
            list.collect { pair -> "${pair[0]}:${pair[1]}" }.join(" ")
        }
        INPUT_QUANTIFICATION_CORRECTED (
            ch_quantification_corrected_unmapped,
            ch_baseline_info,
            //CHECKMATE.out.best_guesses
            ch_checkmate_best_guesses
        )
        ch_quantification_corrected = INPUT_QUANTIFICATION_CORRECTED.out.quantification_corrected_input

        QUANTIFICATION_CORRECTED (
            ch_quantification_corrected
        )
        ch_quantification_corrected_virus_input = QUANTIFICATION_CORRECTED.out.bam_alluniquedcs.join(QUANTIFICATION_CORRECTED.out.bai_alluniquedcs)
        ch_quantification_corrected_virus_input = ch_quantification_corrected_virus_input
        .map { meta, bam, bai -> tuple( meta + [ consensus: "all.unique.dcs"], bam, bai ) }
        ch_genotype_corrected_fasta_fai = QUANTIFICATION_CORRECTED.out.fasta_fai_genotype

        // Report: Quantification coverage (with -f 2) on baseline-corrected genotype
        COVERAGE_QUANTIFICATION_CORRECTED_VIRUS (
            ch_quantification_corrected_virus_input,
            ch_genotype_corrected_fasta_fai,
            QUANTIFICATION_CORRECTED.out.dict_genotype,
            QUANTIFICATION_CORRECTED.out.bed_genotype
        )
        ch_coverage_quantification_corrected_virus = COVERAGE_QUANTIFICATION_CORRECTED_VIRUS.out.summary.collect()
        ch_read_families_corrected_virus = QUANTIFICATION_CORRECTED.out.read_families.collect()

        // Report: fragment insert size virus on corrected genotype
        if ( !params.skip_insertsize ){
            ch_insertsize_corrected = ch_quantification_corrected_virus_input
            .map { meta, bam, bai ->
                tuple( meta + [ insertsize: "corrected" ], bam, bai )
            }
            INSERTSIZE_CORRECTED (
                ch_insertsize_corrected
            )
            ch_insertsize_report_corrected_virus  = INSERTSIZE_CORRECTED.out.metrics.map { meta, file -> file }.collect()
        }
    }

    //
    // MODULE: Run SearcHPV for HPV integration on corrected genotypes
    //

    if (!params.skip_integration && !params.skip_genotype_correction) {
        REFORMAT_FASTQ2 (
            ch_tag2header_reads
        )
        ch_tag2header_reads_v2 = REFORMAT_FASTQ2.out.reads 
        .map { meta, r1, r2 ->
            def cleaned_meta = meta.subMap(['id'])
            tuple( cleaned_meta, [ r1, r2 ] ) 
        }
        ch_quantification_v2 = ch_quantification_corrected
        .map { meta, bam_unmapped, genotype ->
            def cleaned_meta = meta.subMap(['id'])
            tuple( cleaned_meta, genotype ) 
        }
        ch_integration_input = ch_quantification_v2
        .join(ch_tag2header_reads_v2)
        .map { meta, genotype, reads ->
            tuple( meta, reads, genotype ) 
        }
        ch_integration_input.view { it -> 
            println(it)
            println(it*.getClass())
        }
        ch_genotype_bwa_index = ch_integration_input
        .map { meta, reads, genotype ->
            def virus_name = ""
            if (genotype instanceof Path || genotype instanceof File) {
                def lines = genotype.readLines()
                virus_name = lines ? lines[0].trim() : null
            } else {
                virus_name = genotype.toString().trim()
            }
            def index = "${params.ref_path_virus}/${virus_name}/"
            tuple( meta + [ virus: virus_name ], file(index, checkIfExists: true) ) 
        }
        ch_bwa_index.view { it -> 
            println(it)
            println(it*.getClass())
        }
        ch_genotype_bwa_index.view { it -> 
            println(it)
            println(it*.getClass())
        }
        ch_integration_input = ch_integration_input.map { meta, reads, genotype -> tuple( meta, reads ) }
        INTEGRATION_SEARCHPV (
            ch_integration_input,
            ch_bwa_index,
            ch_genotype_bwa_index 
        ) 
        ch_integration_report = INTEGRATION_SEARCHPV.out.HPVfusionPointContigSrNum.collect() 
    }

    //
    // Summary: on-target rate (hsmetrics), coverage qc, genotyping, coverage quantification, saturation rate
    //
    ch_coverage_qc_bwa_report           = COVERAGE_QC_BWA.out.summary.collect()
    ch_coverage_qc_dcs_reprot           = COVERAGE_QC_DCS.out.summary.collect()
    ch_coverage_qc_dcssc_report         = COVERAGE_QC_DCSSC.out.summary.collect()
    ch_coverage_qc_sscs_report          = COVERAGE_QC_SSCS.out.summary.collect()
    ch_coverage_qc_sscssc_report        = COVERAGE_QC_SSCSSC.out.summary.collect()
    ch_coverage_qc_alluniquedcs_report  = COVERAGE_QC_ALLUNIQUEDCS.out.summary.collect()

    genotypes_info       = ch_genotypes.map { map, file -> file}

    ch_coverage_quantification_hg              = COVERAGE_QUANTIFICATION_HG.out.summary.collect()
    ch_coverage_quantification_virus           = COVERAGE_QUANTIFICATION_VIRUS.out.summary.collect()

    ch_read_families                  = CONSENSUSCRUNCHER.out.optional_read_families.collect()
    ch_read_families_virus            = QUANTIFICATION_DOMINANT.out.read_families.collect()

    SUMMARIZE_REPORTS (
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
    if ( !params.skip_genotype_correction ) {
        SUMMARIZE_QUANTIFICATION_CORRECTED (
	    ch_coverage_quantification_hg, 
	    ch_coverage_quantification_virus,
	    ch_read_families, 
	    ch_read_families_corrected_virus
        )     
    }
    if ( !params.skip_insertsize ) {
        SUMMARIZE_INSERTSIZE (
            ch_insertsize_report_hg,
            ch_insertsize_report_dominant_virus,
            ch_insertsize_report_corrected_virus
        )     
    }
    if ( !params.skip_integration ) {
        SUMMARIZE_INTEGRATION (
            ch_integration_report
        )     
    }

    emit:
    cc_dir                  = CONSENSUSCRUNCHER.out.consensus_dir
    cc_dir_genotype         = QUANTIFICATION_DOMINANT.out.consensus_dir
    multiqc_report          = "Later"
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
