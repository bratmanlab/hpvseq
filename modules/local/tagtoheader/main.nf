process TAGTOHEADER {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    path(blist)

    output:
    tuple val(meta), path("*_barcode_R1.fastq.gz"), path("*_barcode_R2.fastq.gz"), emit: reads

    script:
    """
    python3 /cluster/home/jfzou/ConsensusCruncher3/ConsensusCruncher/extract_barcodes_overlapStartswith.py --read1 ${reads[0]} --read2 ${reads[1]} --outfile $meta.id --blist ${blist}
    gzip ${meta.id}_barcode_R1.fastq
    gzip ${meta.id}_barcode_R2.fastq
    """

    stub:
    """
    touch ${meta.id}_barcode_R1.fastq.gz
    touch ${meta.id}_barcode_R2.fastq.gz
    """
}

