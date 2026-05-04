process SUMMARIZE_QUANTIFICATION_CORRECTED {
    tag "summarize_quantification_corrected"
    label 'process_single'

    input:
    path coverage_quantification_hg_files 
    path coverage_quantification_corrected_virus_files
    path read_families_hg_files
    path read_families_virus_corrected_files

    output:
    path("*coverage.quantification.corrected.summ.txt"), emit: summ_coverage_quantification_corrected 
    path("*coverage.quantification.corrected.summ.table.txt"), emit: summ_coverage_quantification_corrected_table
    path("*saturation.corrected.summ.table.txt"), emit: summ_saturation_corrected_table 

    script:
    summ_coverage_quantification_corrected_file = "coverage.quantification.corrected.summ.txt"
    summ_coverage_quantification_corrected_table_file = "coverage.quantification.corrected.summ.table.txt"
    summ_saturation_corrected_table_file = "saturation.corrected.summ.table.txt"
    """
    ## coverage
    get_cov(){
        local infile=\$1
        local label=\$2
        local outfile=\$3

	samp=\$(head -2 \$infile | tail -n 1| awk '{print \$1}')
        bamtype=\$label
	depth=\$(head -n 2 \$infile | tail -n 1 | awk '{print \$3}')
	echo "\$samp \$bamtype \$depth" >> \$outfile
    }
    pivot() {
        ## Sample, Column, Value
        local infile=\$1
        local outfile=\$2
        awk '{print \$1}' \$infile | sort | uniq >| samps
        awk '{print \$2}' \$infile | sort | uniq >| cols
        header="Sample"
        cat cols | while read col;do header=\$header" "\$col; echo \$header;done | tail -1 > \$outfile
    
        last=\$(tail -1 "cols")
        cat "samps" | while read samp;do
            row=\$samp
            cat "cols" | while read col;do
                n=\$(awk -v col=\$col -v samp=\$samp '{if(\$1==samp && \$2==col) print \$3}' \$infile)
                if [ "\$n" == "" ];then
                    n=0
                fi
                row=\$row" "\$n
                if [ "\$col" == "\$last" ];then
                    echo \$row
                fi
                done
            done >> \$outfile
    }

    ## coverage quantification on baseline-corrected genotype
    if [ -f "${summ_coverage_quantification_corrected_file}" ];then rm ${summ_coverage_quantification_corrected_file};fi 
    for file in ${coverage_quantification_hg_files};do get_cov \$file "hg" ${summ_coverage_quantification_corrected_file};done 
    for file in ${coverage_quantification_corrected_virus_files};do get_cov \$file "virus" ${summ_coverage_quantification_corrected_file};done 
    pivot ${summ_coverage_quantification_corrected_file} ${summ_coverage_quantification_corrected_table_file}

    ## saturation rate
    saturation() {
        local infile=\$1
        local label=\$2
        local outfile_total=\$3
        local outfile_nonsingles=\$4
        local outfile_rate=\$5

        samp=\$(echo \$infile | awk -F"_${params.genome}" '{print \$1}')
        total=\$(awk 'NR>1{s+=\$2}END{print s}' \$infile)
        singles=\$(awk '{if(\$1==1) print \$2}' \$infile)
        nonsingles=\$((\$total - \$singles))
        rate=\$(echo "scale=4; \$nonsingles/\$total*100" | bc)
        echo "\$samp \$label \$total" >> \$outfile_total
        echo "\$samp \$label \$nonsingles" >> \$outfile_nonsingles
        echo "\$samp \$label \$rate" >> \$outfile_rate
    }

    for file in tmp_saturation_total tmp_saturation_nonsingles tmp_saturation_rate;do
        if [ -f "\$file" ];then rm \$file;fi
        done

    for file in ${read_families_hg_files};do
        saturation \$file hg tmp_saturation_total tmp_saturation_nonsingles tmp_saturation_rate
        done
    for file in ${read_families_virus_corrected_files};do
        saturation \$file virus tmp_saturation_total tmp_saturation_nonsingles tmp_saturation_rate
        done

    for file in tmp_saturation_total tmp_saturation_nonsingles tmp_saturation_rate;do
        pivot \$file \${file}_table
        done

    for f in tmp_saturation_nonsingles_table tmp_saturation_rate_table; do
        awk 'NR==FNR { a[i++]=\$1; next } { d[\$1]=\$0 } END { for(j=0; j<i; j++) print d[a[j]] }' tmp_saturation_total_table \$f > \${f}_2
    done

    awk -v label="families_Total" 'NR==1{print \$1, label"_"\$2, label"_"\$3}NR>1{print \$0}' tmp_saturation_total_table >| tmp_saturation_total_table_2
    awk -v label="families_gt2" 'NR==1{print label"_"\$2, label"_"\$3}NR>1{print \$0}' tmp_saturation_nonsingles_table_2 >| tmp_saturation_nonsingles_table_3
    awk -v label="PCT_families_gt2" 'NR==1{print label"_"\$2, label"_"\$3}NR>1{print \$0}' tmp_saturation_rate_table_2 >| tmp_saturation_rate_table_3
    paste tmp_saturation_total_table_2 tmp_saturation_nonsingles_table_3 tmp_saturation_rate_table_3 > ${summ_saturation_corrected_table_file}
    #rm tmp_saturation_* 
    """

    stub:
    """
    touch coverage.quantification.corrected.summ.txt
    touch coverage.quantification.corrected.summ.table.txt
    """
}

