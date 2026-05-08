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
    bam           // channel: [ val(meta), path(bam), path(bai) ]
    index         // channel: path(bwa index)
    fasta         // channel: [ val(meta2), path(fasta), path(fai) ]
    cytoband      // channel: [ val(meta2), path(cytoband) ]
    genotypes     // channel: [ val(meta3), path(genotypes) ]

    main:

    //
    // Unmapped reads
    //
    SAMTOOLS_VIEW_UNMAPPED ( 
       bam,
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
    SAMTOOLS_COLLATEFASTQ (
        ch_unmapped,
        fasta.first(),
        false
    )

    //
    // Align to the given genome 
    //
    ch_fastq = SAMTOOLS_COLLATEFASTQ.out.fastq
    BWA_MEM (
        ch_fastq,
        index.map { index_path -> [ [id: params.genome_genotyping], index_path ] },
        [ [:], [] ],    // No fasta needed for BAM output
        false           // sort_bam - we'll sort with samtools for consistency
    )
    ch_orig_bam = BWA_MEM.out.bam
    ch_orig_bam = ch_orig_bam
        .combine(fasta)
        .map { meta, bam, meta_genome, fasta, fai ->
            [ meta + [ ref2: meta_genome.id, type2: "Aligned", consensus: "none"], bam ]
        }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS (
        ch_orig_bam,
        ch_orig_bam
        .combine(fasta)
        .map { meta, bam, meta_genome, fasta, fai ->
            [ meta, fasta, fai]
        }.first()
    )

    //
    // ConsensusCruncher: genotyping
    //
    ch_cc_genotyping = BAM_SORT_STATS_SAMTOOLS.out.bam
       .join(BAM_SORT_STATS_SAMTOOLS.out.index)


    CONSENSUSCRUNCHER (
        ch_cc_genotyping,
        cytoband
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
        [ meta + [ type: "Unmapped", ref2: "virus", type2: "Unknown", consensus: "none" ], file ]
    }

    ch_quantification = ch_unmapped
        .join(ch_genotype, remainder: true)
        .branch { meta, bam_unmapped, genotype_file ->
            //def f = genotype_file instanceof Path ? genotype_file : file(genotype_file)
            has_data: genotype_file != null && genotype_file.exists() && genotype_file.size() > 0
            empty   : true
        }
 

    emit:
    quantification_input     = ch_quantification.has_data
    genotyping_reads   = GENOTYPING_READS.out.genotyping_reads
    bam                = BAM_SORT_STATS_SAMTOOLS.out.bam
    bai                = BAM_SORT_STATS_SAMTOOLS.out.index
    consensus_dir      = CONSENSUSCRUNCHER.out.consensus_dir
    read_families      = CONSENSUSCRUNCHER.out.optional_read_families
    bam_alluniquedcs   = CONSENSUSCRUNCHER.out.bam_alluniquedcs
    bai_alluniquedcs   = CONSENSUSCRUNCHER.out.bai_alluniquedcs
    bam_unmapped       = ch_unmapped
}

