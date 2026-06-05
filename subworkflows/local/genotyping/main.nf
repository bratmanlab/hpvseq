//
// Genotyping
//
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_UNMAPPED  } from '../../../modules/nf-core/samtools/view'
include { BWA_MEM                     } from '../../../modules/nf-core/bwa/mem'
include { CONSENSUSCRUNCHER           } from '../../../modules/local/consensuscruncher'
include { SAMTOOLS_COLLATEFASTQ       } from '../../../modules/nf-core/samtools/collatefastq'  
include { GENOTYPING_READS            } from '../../../modules/local/genotypingreads'  

include { BAM_SORT_STATS_SAMTOOLS     } from '../../nf-core/bam_sort_stats_samtools'

workflow GENOTYPING {
    take:
    input           // channel: [ val(meta), path(bam), path(bai), path(index), path(fasta, path(fai), path(cytoband)]
    //bam           // channel: [ val(meta), path(bam), path(bai)]
    //index         //  value: path(bwa index)
    //fasta         //  value: path(fasta)
    //fai         //    value: path(fai) 
    //cytoband      //  value: path(cytoband)
    genotypes     // channel: [ val(meta2), path(genotypes) ]

    main:

    //
    // Unmapped reads
    //
    ch_bam = input.map { meta, bam, bai, index, fasta, fai, cytoband -> tuple( meta, bam, bai ) }
    SAMTOOLS_VIEW_UNMAPPED ( 
       ch_bam,
       [[:], [], []],
       [[:], []],
       [[:], []],
       ""
    )
    ch_unmapped = SAMTOOLS_VIEW_UNMAPPED.out.bam
    .map { meta, bam -> 
        [ meta + [ type: "Unmapped"], bam ]
    }

    //
    // Bam to fastq 
    //
/*
    ch_fasta_fai = input
    .map { meta, bam, bai, index, fasta, fai, cytoband ->
        tuple( meta, fasta, fai )
    }


    ch_fasta_fai = bam
    .map { meta, bam, bai ->
        tuple( meta, fasta, fai )
    }
*/
    ch_fasta_fai = input
    .map { meta, bam, bai, index, fasta, fai, cytoband -> 
        tuple( meta + [ type: "Unmapped"], fasta, fai ) 
    }
    ch_collatefastq = ch_unmapped.join(ch_fasta_fai)
    //ch_unmapped.view { it -> println("Genotyping"); println(it); println(it*.getClass()) }
    SAMTOOLS_COLLATEFASTQ (
        ch_collatefastq,
        //ch_unmapped,
        //fasta.first(),
        //ch_fasta_fai,
        false
    )

    //
    // Align to the given genome 
    //
    //ch_index = channel.value([ [id: params.genome_genotyping], index ])
    ch_index = input.map { meta, bam, bai, index, fasta, fai, cytoband -> tuple( meta + [ type: "Unmapped"], index ) }
    ch_fastq = SAMTOOLS_COLLATEFASTQ.out.fastq
    .join(ch_index)
    .map { meta, reads, index -> 
        tuple( meta, reads, index, [], [])
    }
    BWA_MEM (
        ch_fastq,
        //index.map { index_path -> [ [id: params.genome_genotyping], index_path ] },
        //ch_index,
        //[ [:], [] ],    // No fasta needed for BAM output
        false           // sort_bam - we'll sort with samtools for consistency
    )
    ch_orig_bam = BWA_MEM.out.bam
    ch_orig_bam = ch_orig_bam
        .map { meta, bam ->
            [ meta + [ ref2: params.genome_genotyping, type2: "Aligned", consensus: "none"], bam ]
        }

    //ch_orig_bam.view { it -> println("Genotyping BWA_MEM"); println(it); println(it*.getClass()) }
    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_fasta_fai = ch_fasta_fai
    .map { meta, fasta, fai ->
        tuple( meta + [ ref2: params.genome_genotyping, type2: "Aligned", consensus: "none"], fasta, fai )
    }
    ch_bam_sort_stats = ch_orig_bam.join(ch_fasta_fai)
    BAM_SORT_STATS_SAMTOOLS (
        ch_bam_sort_stats
    //    ch_orig_bam,
    //    ch_fasta_fai
    )
    //BAM_SORT_STATS_SAMTOOLS.out.bam.view { it -> println("Genotyping BAM_SORT_STATS_SAMTOOLS"); println(it); println(it*.getClass()) }
    
    //
    // ConsensusCruncher: genotyping
    //
    ch_cytoband = input
    .map { meta, bam, bai, index, fasta, fai, cytoband -> 
        tuple( meta + [ type: "Unmapped", ref2: params.genome_genotyping, type2: "Aligned", consensus: "none"], cytoband ) 
    }
    //ch_cytoband.view { it -> println("Genotyping cytoband"); println(it); println(it*.getClass()) }
    ch_cc_genotyping = BAM_SORT_STATS_SAMTOOLS.out.bam
    .join(BAM_SORT_STATS_SAMTOOLS.out.index)
    .join(ch_cytoband)
    .map { meta, bam, bai, cytoband ->
        tuple( meta + [ cc_type: "genotyping" ], bam, bai, cytoband ) 
    }
    //ch_cc_genotyping.view { it -> println("Genotyping ConsensusCruncher"); println(it); println(it*.getClass()) }
    CONSENSUSCRUNCHER (
        ch_cc_genotyping
        //cytoband
    )
    ch_alluniqdcs = CONSENSUSCRUNCHER.out.bam_alluniquedcs
        .join(CONSENSUSCRUNCHER.out.bai_alluniquedcs)

    //
    // Output genotypes and put dominant as the first 
    //
    GENOTYPING_READS (
       ch_alluniqdcs,
       params.mapq,
       genotypes 
    )
 
    ch_genotype = GENOTYPING_READS.out.genotype
    .map { meta, file ->
        def new_meta = meta.clone()
        new_meta.remove("cc_type")
        tuple( [ new_meta + [ type: "Unmapped", ref2: "virus", type2: "Unknown", consensus: "none" ], file ] )
    }
    ch_quantification = ch_unmapped
        .join(ch_genotype, remainder: true)
        .branch { meta, bam_unmapped, genotype_file ->
            //def f = genotype_file instanceof Path ? genotype_file : file(genotype_file)
            //has_data: genotype_file != null && genotype_file.exists() && genotype_file.size() > 0
            has_data: genotype_file instanceof Path && genotype_file.exists() && genotype_file.size() > 0
            empty   : true
        }
    ch_quantification_input = ch_quantification.has_data
    .map { meta, bam_unmapped, genotype_file ->
        def lines = genotype_file.readLines()
        virus_name = lines ? lines[0].trim() : null
        tuple( meta + [ virus: virus_name ], bam_unmapped )
    } 

    emit:
    quantification_input     = ch_quantification_input
    genotyping_reads   = GENOTYPING_READS.out.genotyping_reads
    bam                = BAM_SORT_STATS_SAMTOOLS.out.bam
    bai                = BAM_SORT_STATS_SAMTOOLS.out.index
    consensus_dir      = CONSENSUSCRUNCHER.out.consensus_dir
    read_families      = CONSENSUSCRUNCHER.out.optional_read_families
    bam_alluniquedcs   = CONSENSUSCRUNCHER.out.bam_alluniquedcs
    bai_alluniquedcs   = CONSENSUSCRUNCHER.out.bai_alluniquedcs
    bam_unmapped       = ch_unmapped
}

