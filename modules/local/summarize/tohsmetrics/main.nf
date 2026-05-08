process SUMMARIZE_TOHSMETRICS {
    tag "summarize_toHsmetrics"
    label 'process_single'

    input:
    path nreads_raw_files
    path nreads_goodumi_files
    path hsmetrics_files

    output:
    path("*nreads.summ.txt"), emit: summ_nreads 
    path("*nreads.summ.table.txt"), emit: summ_nreads_table 
    path("*hsmetrics.summ.txt"), emit: summ_hsmetrics 
    path("*hsmetrics.summ.table.txt"), emit: summ_hsmetrics_table
    path("summ_tohsmetrics.txt"), emit: summ_tohsmetrics 

    script:
    hsmetrics_file1 = hsmetrics_files[0]
    summ_nreads_file = "nreads.summ.txt"
    summ_nreads_table_file = "nreads.summ.table.txt"
    summ_hsmetrics_file = "hsmetrics.summ.txt"
    summ_hsmetrics_table_file = "hsmetrics.summ.table.txt"
    """
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

    ## nreads: raw, goodumi
    for file in ${nreads_raw_files};do
        samp=\${file/_raw*}
        n=\$(head -1 \$file)
        echo "\$samp raw \$n"
        done > ${summ_nreads_file} 
    for file in ${nreads_goodumi_files};do
        samp=\${file/_goodUMI*}
        n=\$(head -1 \$file)
        echo "\$samp goodUMI \$n"
        done >> ${summ_nreads_file} 
         
    pivot ${summ_nreads_file} ${summ_nreads_table_file}

    ## on-target rate
    head ${hsmetrics_file1} -n 7 | tail -n 1 | awk -F'\t' '{print "samp",\$2,\$21,\$3,\$23,\$31,\$32,\$28,\$29,\$4,\$5,\$6,\$30,\$8,\$10,\$34,\$11,\$12}' > ${summ_hsmetrics_file}
    for file in ${hsmetrics_files};do
        samp=\$(echo \$file | awk -F'.CollectHsMetrics' '{print \$1}')
        head \$file -n 8 | tail -n 1 | awk -F'\t' -v samp=\$samp '{print samp,\$2,\$21,\$3,\$23,\$31,\$32,\$28,\$29,\$4,\$5,\$6,\$30,\$8,\$10,\$34,\$11,\$12}' >> ${summ_hsmetrics_file}
        done

    awk 'NR==1{print "Sample", "PCT_OFF_BAIT", "PCT_NEAR_BAIT", "PCT_ON_BAIT", "PCT_ON_TARGET"}NR>1{offbait=\$12/\$8*100; nearbait=\$11/\$8*100; onbait=\$10/\$8*100; ontarget=\$13/\$8*100; print \$1, offbait, nearbait, onbait, ontarget}' ${summ_hsmetrics_file} > ${summ_hsmetrics_table_file}

    ## merge nreads and on-target rate
    for f in ${summ_hsmetrics_table_file}; do
        awk 'NR==FNR { a[i++]=\$1; next } { d[\$1]=\$0 } END { for(j=0; j<i; j++) print d[a[j]] }' ${summ_nreads_table_file} \$f > \${f}_2
    done

    paste -d " " ${summ_nreads_table_file} ${summ_hsmetrics_table_file}_2 > summ_tohsmetrics.txt

    """

    stub:
    """
    touch nreads.summ.txt
    touch nreads.summ.table.txt
    touch hsmetrics.summ.txt
    touch hsmetrics.summ.table.txt
    """
}

