process INTEGRATION_SEARCHPV {
    tag "$meta.id"
    label 'process_higher'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index_hg)
    tuple val(meta3), path(index_virus)

    output:
    tuple val(meta), path("${meta.id}__${meta3.virus}"), emit: searchpv_dir
    path("${meta.id}__${meta3.virus}/call_fusion_virus/${meta.id}_HPVfusionPointContig.txt"), optional: true, emit: HPVfusionPointContigSrNum

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-mem ${task.memory.giga} -nt ${task.cpus}'
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    echo \$INDEX
    fasta_virus="${meta3.virus}.fasta"
    for index in \${INDEX[@]};do
        if [[ "\$index" == *"\$fasta_virus" ]]; then
            ref_virus=\$index
        fi    
        if [[ "\$index" != *"\$fasta_virus" ]]; then
            ref_hg=\$index
        fi    
        done

    echo \$ref_hg
    echo \$ref_virus
    searcHPV -fastq1 ${reads[0]} \
        -fastq2 ${reads[1]} \
	-humRef \${ref_hg} \
	-virRef \${ref_virus} \
	-output ${meta.id}__${meta3.virus} \
        ${args}

    ## Collect srNum
    txt="${meta.id}__${meta3.virus}/call_fusion_virus/HPVfusionPointContig.txt"
    n=\$(wc -l \$txt | awk '{print \$1}')
    if [[ "\$n" -gt "1" ]];then
        Dir=\$(dirname \$txt);
        awk '{if(\$0!="") print \$0}' \$txt > tmp1_${meta.id}
        echo "srNum" > tmp2_${meta.id}
        awk 'NR>1{print \$1}' tmp1_${meta.id} | while read contig;do
            site=\${contig/.Contig*};
            sr=\$(awk -v contig=\$contig '{if(\$1==contig) print \$2}' \$Dir/\$site/srNum.txt);
            echo \$sr
            done >> tmp2_${meta.id}
        paste tmp1_${meta.id} tmp2_${meta.id} > "${meta.id}__${meta3.virus}/call_fusion_virus/${meta.id}_HPVfusionPointContig.txt"
        #rm tmp*_${meta.id}
    fi
    """

    stub:
    """
    mkdir -p ${meta.id}__${meta3.virus}/call_fusion_virus
    touch ${meta.id}__${meta3.virus}/call_fusion_virus/HPVfusionPointContig.txt
    """
}

