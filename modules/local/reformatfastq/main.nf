process REFORMAT_FASTQ {
    tag "${meta.id}_REFORMATFASTQ"
    label 'process_single'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.reformatted_1.fastq.gz"), path("*.reformatted_2.fastq.gz"), emit: reads

    // >=2 fields but first two are correct
    // zcat ${reads[0]} | awk '{if(NR%4==1) print \$1, \$2; else if(NR%4==3) print "+"; else print \$0}' | gzip > ${meta.id}.reformatted_1.fastq.gz
    // zcat ${reads[1]} | awk '{if(NR%4==1) print \$1, \$2; else if(NR%4==3) print "+"; else print \$0}' | gzip > ${meta.id}.reformatted_2.fastq.gz
    // >=2 fields but second includes flowcell information
    script:
    """
    zcat ${reads[0]} | awk '{if(NR%4==1) {split(\$1, arr1, ":"); split(\$2, arr2, ":"); if(length(arr1)>1) print \$1, \$2; else if (length(arr2)>1) print "@"\$2, \$3; else print \$1, \$2} else if(NR%4==3) print "+"; else print \$0}' | gzip > ${meta.id}.reformatted_1.fastq.gz
    zcat ${reads[1]} | awk '{if(NR%4==1) {split(\$1, arr1, ":"); split(\$2, arr2, ":"); if(length(arr1)>1) print \$1, \$2; else if (length(arr2)>1) print "@"\$2, \$3; else print \$1, \$2} else if(NR%4==3) print "+"; else print \$0}' | gzip > ${meta.id}.reformatted_2.fastq.gz
    """

    stub:
    """
    touch ${meta.id}.reformatted_1.fastq.gz
    touch ${meta.id}.reformatted_2.fastq.gz
    """
}

