process CORRECT_GENOTYPE {
    label 'process_single'

    input:
    val(baseline_info)
    path(mates)
    val(others_info)

    output:
    path("corrected_genotypes.txt"), emit: corrected_genotype
    path("mates.txt"), emit: mates_clean 

    script:
    """
    # mates
    n=\$(wc -l ${mates} | awk '{print \$1}')
    for i in \$(seq 1 \$n);do 
        awk -F'\\t' -v line=\$i '
        NR==line{
            gsub(/\\\133/, ""); 
            gsub(/\\047/, ""); 
            s=\$1; 
            for(j=2;j<=NF;j++){
                split(\$j, arr, ","); 
                s=s" "arr[1]
            }
            print s
        }
        ' ${mates}
        done > mates.txt
    for line in ${baseline_info};do 
        samp=\$(echo \$line | cut -d ":" -f1)
        genotype=\$(echo \$line | cut -d ":" -f2)
        awk -v samp=\$samp '{if(\$1==samp".bam") print \$0}' mates.txt | awk -v genotype=\$genotype '{gsub(/.bam\$/, ""); gsub(/.bam /, " "); for(i=1;i<=NF;i++){print \$i, genotype }}' 
        done > corrected_genotypes.txt 
    for line in ${others_info};do 
        samp=\$(echo \$line | cut -d ":" -f1)
        genotype=\$(echo \$line | cut -d ":" -f2)
        n=\$(awk -v samp=\$samp 'BEGIN{n=0}{if(\$1==samp) n+=1}END{print n}' corrected_genotypes.txt )
        if [[ "\$n" == "0" ]];then
            echo \$samp \$genotype
        fi
        done >> corrected_genotypes.txt 
    """

    stub:
    """
    touch corrected_genotypes.txt
    """
}

