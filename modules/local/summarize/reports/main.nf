process SUMMARIZE_REPORTS {
    tag "summarize_report"

    input:
    path coverage_bwa_files
    path coverage_dcs_files
    path coverage_dcssc_files
    path coverage_sscs_files
    path coverage_sscssc_files
    path coverage_alluniquedcs_files
    path genotyping_files
    path genotypes
    path coverage_quantification_hg_files 
    path coverage_quantification_virus_files 
    path read_families_hg_files 
    path read_families_virus_files 

    output:
    path("*coverage.qc.summ.txt"), emit: summ_coverage_qc 
    path("*coverage.qc.summ.table.txt"), emit: summ_coverage_qc_table 
    path("*genotyping.summ.table.txt"), optional: true, emit: summ_genotyping_table 
    path("*coverage.quantification.summ.txt"), emit: summ_coverage_quantification 
    path("*coverage.quantification.summ.table.txt"), emit: summ_coverage_quantification_table 
    path("*saturation.summ.table.txt"), emit: summ_saturation_table 

    script:
    coverage_bwa_file1 = coverage_bwa_files[0]
    genotyping_file1 = genotyping_files[0]
    summ_hsmetrics_file = "hsmetrics.summ.txt"
    summ_hsmetrics_table_file = "hsmetrics.summ.table.txt"
    summ_coverage_qc_file = "coverage.qc.summ.txt"
    summ_coverage_qc_table_file = "coverage.qc.summ.table.txt"
    summ_genotyping_file = "genotyping.summ.txt"
    summ_genotyping_table_file = "genotyping.summ.table.txt"
    summ_coverage_quantification_file = "coverage.quantification.summ.txt"
    summ_coverage_quantification_table_file = "coverage.quantification.summ.table.txt"
    summ_saturation_table_file = "saturation.summ.table.txt"
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

    ##echo "Sample Type Depth" > ${summ_coverage_qc_file}
    if [ -f "${summ_coverage_qc_file}" ];then rm ${summ_coverage_qc_file};fi 
    for file in ${coverage_bwa_files};do get_cov \$file "depth_qc_bwa" ${summ_coverage_qc_file};done 
    for file in ${coverage_dcs_files};do get_cov \$file "depth_qc_dcs" ${summ_coverage_qc_file};done 
    for file in ${coverage_dcssc_files};do get_cov \$file "depth_qc_dcs.sc" ${summ_coverage_qc_file};done 
    for file in ${coverage_sscs_files};do get_cov \$file "depth_qc_sscs" ${summ_coverage_qc_file};done
    for file in ${coverage_sscssc_files};do get_cov \$file "depth_qc_sscs.sc" ${summ_coverage_qc_file};done 
    for file in ${coverage_alluniquedcs_files};do get_cov \$file "depth_qc_all.unique.dcs" ${summ_coverage_qc_file};done 
    pivot ${summ_coverage_qc_file} ${summ_coverage_qc_table_file}

    ## genotyping
    for file in ${genotyping_files};do
        samp=\$(echo \$file | awk -F'_${params.genome_genotyping}' '{print \$1}')
        awk -v samp=\${samp} '{split(\$1, arr, "|"); print samp, arr[1]"|"arr[2], \$2}' \$file >> ${summ_genotyping_file}
        done

    if [ -f "${summ_genotyping_file}" ] && [ -s "${summ_genotyping_file}" ];then
        awk 'NR==FNR{pat[\$1]=\$2; next} \$2 in pat{print \$1, pat[\$2], \$3}' ${genotypes} ${summ_genotyping_file} > tmp_${summ_genotyping_file}
        pivot tmp_${summ_genotyping_file} ${summ_genotyping_table_file}
    fi

    ## coverage quantification
    if [ -f "${summ_coverage_quantification_file}" ];then rm ${summ_coverage_quantification_file};fi 
    for file in ${coverage_quantification_hg_files};do get_cov \$file "hg" ${summ_coverage_quantification_file};done 
    for file in ${coverage_quantification_virus_files};do get_cov \$file "virus" ${summ_coverage_quantification_file};done 
    pivot ${summ_coverage_quantification_file} ${summ_coverage_quantification_table_file}

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
    for file in ${read_families_virus_files};do
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

    paste -d " " tmp_saturation_total_table_2 tmp_saturation_nonsingles_table_3 tmp_saturation_rate_table_3 > ${summ_saturation_table_file}
    #rm tmp_saturation_* 
    """

    stub:
    """
    touch coverage.qc.summ.txt
    touch coverage.qc.summ.table.txt
    """
}

