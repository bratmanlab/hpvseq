process REFORMAT_FASTQ2 {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.reformatted2_1.fastq.gz"), path("*.reformatted2_2.fastq.gz"), emit: reads

    script:
    """
    sample_index="GCTATCCT+TGCGTAGA"
    zcat ${reads[0]} | awk -v sample_index=\$sample_index '{if(NR%4==1 && \$1~/^@/){split(\$1, arr, "|"); print arr[1], "1:N:0:"sample_index}else{print \$0}}' | gzip > ${meta.id}.reformatted2_1.fastq.gz 
    zcat ${reads[1]} | awk -v sample_index=\$sample_index '{if(NR%4==1 && \$1~/^@/){split(\$1, arr, "|"); print arr[1], "2:N:0:"sample_index}else{print \$0}}' | gzip > ${meta.id}.reformatted2_2.fastq.gz 
    """

    stub:
    """
    touch ${meta.id}.reformatted2_1.fastq.gz
    touch ${meta.id}.reformatted2_2.fastq.gz
    """
}

