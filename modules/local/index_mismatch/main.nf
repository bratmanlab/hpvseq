process INDEX_MISMATCH {
    tag "${meta.id}_INDEXMISMATCH"
    label 'process_single'

    input:
    tuple val(meta), path(reads)
    path(index)
    val(index_mismatch)

    output:
    tuple val(meta), path(reads), path("pass.txt"), path("*sample_index.txt"), emit: reads

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut=${index_mismatch}
    samp_index=${meta.id}.R1.index
    zcat ${reads[0]} | awk '\$1~/^@/ && NR%4==1{print \$0}' | awk '{split(\$2, arr, ":"); split(arr[4], arr1, "+"); print arr1[1], arr1[2]}' | sort | uniq -c | awk '{print \$2, \$3, \$1}' | sort -k3,3nr | while read i7_i5;do
	i7=\$(echo \$i7_i5 | awk '{print \$1}')
	i5=\$(echo \$i7_i5 | awk '{print \$2}')
	i7_len=\$(echo \$i7 | awk '{print length(\$1)}')
	i5_len=\$(echo \$i5 | awk '{print length(\$1)}')

        count=0
        while read line;do
            true_i7=\$(echo \$line | awk '{print \$1}')
            true_i5=\$(echo \$line | awk '{print \$2}')
            mismatch_i7=\$(echo \$i7 | awk -v i7_len=\${i7_len} -v true_i7=\${true_i7} '{split(\$0, arr1, ""); split(true_i7, arr2, ""); for(i=1;i<=length(arr1);i++) if(arr1[i]==arr2[i]) n+=1}END{print i7_len-n}')
            mismatch_i5=\$(echo \$i5 | awk -v i5_len=\${i5_len} -v true_i5=\${true_i5} '{split(\$0, arr1, ""); split(true_i5, arr2, ""); for(i=1;i<=length(arr1);i++) if(arr1[i]==arr2[i]) n+=1}END{print i5_len-n}')
            if [[ "\$mismatch_i7" -le "\$cut" ]] && [[ "\$mismatch_i5" -le "\$cut" ]];then
                count=\$(( \$count + 1 ))
            fi
            done < <(sort $index | uniq)
        echo "\$i7_i5 \$count" 
    done > \$samp_index
    
    if [ -f "\$samp_index" ];then
        n=\$(wc -l \$samp_index | awk '{print \$1}')
        npass=\$(awk 'BEGIN{s=0}{if(\$3>0) s+=1}END{print s}' \$samp_index)
        if [[ "\$n" == "\$npass" ]];then
            echo "1" > pass.txt # pass
            head -1 \$samp_index | awk -v samp=${meta.id} '{print \$1, \$2, samp}' > ${prefix}_sample_index.txt
        else
            echo "0" > pass.txt
            echo "null" > ${prefix}_sample_index.txt 
        fi
    else
        echo "0" > pass.txt
        echo "null" > ${prefix}_sample_index.txt 
    fi

    """

    stub:
    """
    touch pass.txt
    """
}

