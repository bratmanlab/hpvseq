//
// Quantification
//
include { BWA_MEM                     } from '../../../modules/nf-core/bwa/mem'
include { CONSENSUSCRUNCHER           } from '../../../modules/local/consensuscruncher'
include { SAMTOOLS_COLLATEFASTQ       } from '../../../modules/nf-core/samtools/collatefastq/main'  

include { BAM_SORT_STATS_SAMTOOLS     } from '../../nf-core/bam_sort_stats_samtools'

workflow QUANTIFICATION {
    take:
    input           // channel: [ val(meta), path(bam), path(genotype) ]

    main:
    bam      = input.map { meta, bam, genotype -> tuple( meta, bam ) } 
    genotype = input.map { meta, bam, genotype -> tuple( meta, genotype ) } 
    //bam.view { meta, bam -> "UNMAPPED BAM PATH : ${bam}" }

    //
    // Prepare virus genomes
    //
/* 
    ch_virus = input
    .map { meta, bam, genotype_file ->
	// 1. Extract the virus name safely from the first line of the file
	// readLines() is okay for tiny files (1-2 lines), but use with caution
	def lines = genotype_file.readLines()
	def virus_name = lines ? lines[0].trim() : null
	tuple( meta, virus_name ) 
    }
*/
    ch_virus = input
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
        tuple( meta, virus_name ) 
    }

    ch_sort_fasta_fai = ch_virus
    .filter { meta, virus_name -> virus_name != null } // Skip if no virus name found
    .map { meta, virus_name ->
	// 2. Construct paths using the extracted name
	def fasta = file("${params.ref_path_virus}/${virus_name}/${virus_name}.fasta")
	def fai   = file("${params.ref_path_virus}/${virus_name}/${virus_name}.fasta.fai")
	tuple( meta + [ virus: virus_name ], fasta, fai ) 
    }

    ch_dict = ch_virus
    .filter { meta, virus_name -> virus_name != null } // Skip if no virus name found
    .map { meta, virus_name ->
	def dict   = file("${params.ref_path_virus}/${virus_name}/${virus_name}.dict")
	tuple( meta, dict ) 
    }
    
    ch_cytoband = ch_virus
    .filter { meta, virus_name -> virus_name != null } // Skip if no virus name found
    .map { meta, virus_name ->
	def cytoband = file("${params.ref_path_virus}/${virus_name}/${virus_name}.cytoband")
	tuple( [id: virus_name], cytoband ) 
    }

    // "" represents whole genome 
    ch_bed = ch_virus
    .map { meta, virus_name ->
	def bedfile   = file("${params.ref_path_virus}/${virus_name}/${virus_name}_E6E7.bed")
        tuple( meta, bedfile.exists() ? bedfile : [] ) 
    }
/*
    ch_cytoband.view { meta, cytoband -> 
        """
        SAMPLE        : ${meta.id}
        CYTOBAND PATH : ${cytoband}
        """
    }
*/
    ch_index = ch_virus
    .filter { meta, virus_name -> virus_name != null } // Skip if no virus name found
    .map { meta, virus_name ->
	def index = file("${params.ref_path_virus}/${virus_name}/")
	tuple( [id: virus_name], index )
    }
/*
    ch_index.view { meta, index -> 
        """
        SAMPLE        : ${meta.id}
        INDEX PATH    : ${index}
        """
    }
*/
    //
    // Bam to fastq 
    //
    SAMTOOLS_COLLATEFASTQ (
        bam,
        [ [:], [], [] ],
        false
    )

    //
    // Align to the given genome 
    //
    ch_fastq = SAMTOOLS_COLLATEFASTQ.out.fastq
/*
    ch_fastq.view { meta, reads -> 
        """
        SAMPLE              : ${meta.id}
        UNMAPPED2FASTQ PATH : ${reads}
        """
    }
*/
    BWA_MEM (
        ch_fastq,
        ch_index,
        [ [:], [] ],    // No fasta needed for BAM output
        false           // sort_bam - we'll sort with samtools for consistency
    )
    ch_orig_bam = BWA_MEM.out.bam
    ch_orig_bam = ch_orig_bam
        .join(ch_virus)
        .map { meta, bam, virus_name ->
            tuple( meta + [ref2: virus_name, type2: "Aligned", consensus: "none"], bam ) 
        }
/*
    ch_orig_bam = ch_orig_bam
        .combine(ch_index)
        .map { meta, bam, meta_genome, index ->
            tuple( meta + [ref2: meta_genome.id, type2: "Aligned", consensus: "none"], bam ) 
        }

    ch_orig_bam.view { meta, bam -> 
        """
        SAMPLE            : ${meta.id}
        BAM PATH RE-ALIGN : ${bam}
        """
    }
*/
    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS (
        ch_orig_bam,
        ch_sort_fasta_fai
    )

    //
    // ConsensusCruncher: genotype
    //
    ch_cc_genotype = BAM_SORT_STATS_SAMTOOLS.out.bam
       .join(BAM_SORT_STATS_SAMTOOLS.out.index)

    CONSENSUSCRUNCHER (
        ch_cc_genotype,
        ch_cytoband
    )

    emit:
    bam                = BAM_SORT_STATS_SAMTOOLS.out.bam
    bai                = BAM_SORT_STATS_SAMTOOLS.out.index
    consensus_dir      = CONSENSUSCRUNCHER.out.consensus_dir
    read_families      = CONSENSUSCRUNCHER.out.optional_read_families
    bam_dcssc          = CONSENSUSCRUNCHER.out.bam_dcssc
    bai_dcssc          = CONSENSUSCRUNCHER.out.bai_dcssc
    bam_alluniquedcs   = CONSENSUSCRUNCHER.out.bam_alluniquedcs
    bai_alluniquedcs   = CONSENSUSCRUNCHER.out.bai_alluniquedcs
    bam_dcs            = CONSENSUSCRUNCHER.out.bam_dcs
    bai_dcs            = CONSENSUSCRUNCHER.out.bai_dcs
    bam_sscs           = CONSENSUSCRUNCHER.out.bam_sscs
    bai_sscs           = CONSENSUSCRUNCHER.out.bai_sscs
    bam_sscssc         = CONSENSUSCRUNCHER.out.bam_sscssc
    bai_sscssc         = CONSENSUSCRUNCHER.out.bai_sscssc
    fasta_fai_genotype = ch_sort_fasta_fai                 //channel
    dict_genotype      = ch_dict                           //channel
    bed_genotype       = ch_bed                            //channel
}

