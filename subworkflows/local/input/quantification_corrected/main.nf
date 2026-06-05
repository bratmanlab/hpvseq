//
// Prepare input for quantification based on baseline-corrected genotype 
//

include { CORRECT_GENOTYPE         } from '../../../../modules/local/correct_genotype'

workflow INPUT_QUANTIFICATION_CORRECTED {
    take:
    bam_unmapped               //      channel: [ val(meta), path(bam) ]
    baseline_ids_genotypes     //          val: "id:genotype id:genotype" for baseline
    mates                      //      channel: path(best_guesses)
    others_ids_genotypes       //          val: "id:genotype id:genotype" for non-baseline

    main:

    //
    // Get corrected genotype 
    //
    CORRECT_GENOTYPE (
        baseline_ids_genotypes,
        mates,
        others_ids_genotypes
    )
    ch_genotypes = CORRECT_GENOTYPE.out.corrected_genotype

    // Generate channel
    ch_corrected_genotypes = ch_genotypes
    .splitText()                 // Break file into individual lines
    .map { line -> 
        def parts = line.trim().split()
        tuple( parts[0], parts[1] ) // This creates the [id, genotype] tuple
    }
    .view { sample, genotype -> "CORRECTED GENOTYPE : $sample -> $genotype" }

    ch_unmapped = bam_unmapped.map { meta, bam ->
        tuple( meta.id, meta, bam ) 
    }
    ch_joined                         = ch_unmapped.join(ch_corrected_genotypes)
    ch_quantification_corrected_input = ch_joined.map { sample, meta, bam, genotype ->
        tuple( meta + [ virus: genotype ], bam)
    }

    emit:
    quantification_corrected_input      = ch_quantification_corrected_input     // channel
    genotypes_corrected                 = ch_genotypes
}
