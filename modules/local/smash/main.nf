process SMASH {

    input:
    path(bams)

    output:
    path("best_guesses.pval_out.txt"), emit: best_guesses

    script:
    """
    ${params.smash} -bam ALL
    """

    stub:
    """
    touch best_guesses.pval_out.txt 
    """
}

