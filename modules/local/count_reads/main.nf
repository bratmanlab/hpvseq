process COUNT_READS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*nreads.txt"), emit: nreads

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat ${reads[0]} | awk '\$1~/^@/ && NR%4==1 {n += 1}END{print n}' > ${prefix}_nreads.txt
    """

    stub:
    """
    touch ${prefix}_nreads.txt 
    """
}

