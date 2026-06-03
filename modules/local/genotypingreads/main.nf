process GENOTYPING_READS {
    tag "${meta.id}_GENOTYPING"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai)
    val(mapq)
    tuple val(meta2), path(genotypes)

    output:
    tuple val(meta), path("*.genotyping.txt"), emit: genotyping_reads, optional: true
    tuple val(meta), path("*.genotype.txt"), emit: genotype, optional: true

    script:
    def file_reads      = "${meta.id}_${params.genome_genotyping}_f2_mapq${mapq}.reads.txt"
    def file_genotyping = "${meta.id}_${params.genome_genotyping}_f2_mapq${mapq}.genotyping.txt"
    """
    # echo "Sample subtype pos" > ${file_reads}
    ${params.samtools_path} view -f 2 -q $mapq $bam > ${file_reads}
    n=\$(wc -l ${file_reads} | awk '{print \$1}')
    if [[ "\$n" -gt "0" ]];then
        awk '{print \$3}' ${file_reads} | sort | uniq -c | sort -k1,1nr | awk '{print \$2, \$1}' > ${file_genotyping}
        genome=\$(head -1 ${file_genotyping} | awk '{split(\$1, arr, "|"); print arr[1]"|"arr[2]}')
        genotype=\$(awk -v genome=\$genome '{if(\$1==genome) print \$2}' $genotypes)
        echo \$genotype > "${meta.id}.genotype.txt"
    fi 
    """

    stub:
    def file_reads      = "${meta.id}_${params.genome_genotyping}_f2_mapq${mapq}.reads.txt"
    def file_genotyping = "${meta.id}_${params.genome_genotyping}_f2_mapq${mapq}.genotyping.txt"
    """
    touch ${file_reads}
    touch ${file_genotyping}
    """
}

