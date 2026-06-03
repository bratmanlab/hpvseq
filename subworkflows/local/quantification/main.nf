//
// Quantification
//
include { BWA_MEM                     } from '../../../modules/nf-core/bwa/mem'
include { CONSENSUSCRUNCHER           } from '../../../modules/local/consensuscruncher'
include { SAMTOOLS_COLLATEFASTQ       } from '../../../modules/nf-core/samtools/collatefastq/main'  

include { BAM_SORT_STATS_SAMTOOLS     } from '../../nf-core/bam_sort_stats_samtools'

workflow QUANTIFICATION {
    take:
    input           // channel: [ val(meta), path(bam) ]

    main:

    //
    // Prepare virus genomes
    //
/*
    ch_bam = input
    .map { meta, bam, genotype ->
        def virus_name = ""
        // Check if genotype is a Path/File or just a String
        if (genotype instanceof Path || genotype instanceof File) {
            // Read first line if it's a file
            def lines = genotype.readLines()
            virus_name = lines ? lines[0].trim() : null
        } else {
            // It's already a string, use it directly
            virus_name = genotype.toString().trim()
        }
        tuple( meta + [ virus: virus_name ], bam ) 
    }
    ch_dict = input
    .map { meta, bam ->
        def virus_name = meta.virus
	def dict   = file("${params.ref_path_virus}/${virus_name}/${virus_name}.dict")
	tuple( meta, dict ) 
    }
*/ 
    ch_cytoband = input 
    .map { meta, bam ->
        def virus_name = meta.virus
	def cytoband = file("${params.ref_path_virus}/${virus_name}/${virus_name}.cytoband")
	tuple( meta, cytoband) 
    }

    //
    // Bam to fastq 
    //
    ch_collatefastq = input.map { meta, bam -> tuple( meta, bam, [], [] )}
    SAMTOOLS_COLLATEFASTQ (
        ch_collatefastq,
        //input,
        //tuple( [:], null, null ),
        //[ [:], [], [] ],
        false
    )

    //
    // Align to the given genome 
    //
    ch_fastq = SAMTOOLS_COLLATEFASTQ.out.fastq
    .map { meta, reads ->
        def virus_name = meta.virus
	def index = file("${params.ref_path_virus}/${virus_name}/")
	tuple( meta, reads, index, [], [] )
    }
    BWA_MEM (
        ch_fastq,
//        ch_index,
//        [ [:], [] ],    // No fasta needed for BAM output
        false           // sort_bam - we'll sort with samtools for consistency
    )
    ch_orig_bam = BWA_MEM.out.bam
    .map { meta, bam ->
        tuple( meta + [ref2: meta.virus, type2: "Aligned", consensus: "none"], bam ) 
    }
    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_bam_sort_stats = ch_orig_bam
    .map { meta, bam ->
	// Construct paths using the extracted name
        def virus_name = meta.virus
	def fasta = file("${params.ref_path_virus}/${virus_name}/${virus_name}.fasta")
	def fai   = file("${params.ref_path_virus}/${virus_name}/${virus_name}.fasta.fai")
	tuple( meta, bam, fasta, fai ) 
    }
    BAM_SORT_STATS_SAMTOOLS (
        ch_bam_sort_stats
        //ch_orig_bam,
        //ch_fasta_fai
    )

    //
    // ConsensusCruncher: genotype
    //
    ch_cc_genotype = BAM_SORT_STATS_SAMTOOLS.out.bam
    .join(BAM_SORT_STATS_SAMTOOLS.out.index)
    .map { meta, bam, bai -> 
        def virus_name = meta.virus
	def cytoband = file("${params.ref_path_virus}/${virus_name}/${virus_name}.cytoband")
	tuple( meta, bam, bai, cytoband) 
    }
    CONSENSUSCRUNCHER (
        ch_cc_genotype
    )
    
    emit:
    bam                = BAM_SORT_STATS_SAMTOOLS.out.bam                // channel: [ meta, bam]
    bai                = BAM_SORT_STATS_SAMTOOLS.out.index              // channel: [ meta, bai]
    consensus_dir      = CONSENSUSCRUNCHER.out.consensus_dir            //     dir:
    read_families      = CONSENSUSCRUNCHER.out.optional_read_families   // channel: [ meta,
    bam_dcssc          = CONSENSUSCRUNCHER.out.bam_dcssc                // channel: [ meta,
    bai_dcssc          = CONSENSUSCRUNCHER.out.bai_dcssc
    bam_alluniquedcs   = CONSENSUSCRUNCHER.out.bam_alluniquedcs
    bai_alluniquedcs   = CONSENSUSCRUNCHER.out.bai_alluniquedcs
    bam_dcs            = CONSENSUSCRUNCHER.out.bam_dcs
    bai_dcs            = CONSENSUSCRUNCHER.out.bai_dcs
    bam_sscs           = CONSENSUSCRUNCHER.out.bam_sscs
    bai_sscs           = CONSENSUSCRUNCHER.out.bai_sscs
    bam_sscssc         = CONSENSUSCRUNCHER.out.bam_sscssc
    bai_sscssc         = CONSENSUSCRUNCHER.out.bai_sscssc
//    fasta_fai_genotype = ch_fasta_fai                                   //channel: [ meta, fasta, fai ]
//    dict_genotype      = ch_dict                                        //channel: [ meta, dict ]
//    bed_genotype       = ch_bed                                         //channel: [ meta, bed ]
}

