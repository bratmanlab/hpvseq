process SUMMARIZE_INSERTSIZE {
    tag "summarize_insertsize"

    input:
    path insertsize_hg_files 
    path insertsize_dominant_virus_files 
    path insertsize_corrected_virus_files 

    output:
    path("*insertsize.hg.summ.txt"), emit: summ_insertsize_hg
    path("*insertsize.dominant.virus.summ.txt"), emit: summ_insertsize_dominant_virus
    path("*insertsize.corrected.virus.summ.txt"), emit: summ_insertsize_corrected_virus 

    script:
    summ_insertsize_hg_file              = "insertsize.hg.summ.txt"
    summ_insertsize_dominant_virus_file  = "insertsize.dominant.virus.summ.txt"
    summ_insertsize_corrected_virus_file = "insertsize.corrected.virus.summ.txt"
    """
    get_data(){
        local infile=\$1
        local outfile=\$2

        samp=\${infile/_${params.genome}*}
        awk -v samp=\$samp 'NR>11{OFS="\t"; print samp, \$1, \$2}' \${infile} >> \$outfile
    }
    ## insert size hg
    echo -e "Sample\tinsert_size\tAll_Reads.fr_count" > ${summ_insertsize_hg_file}
    for file in ${insertsize_hg_files};do
        get_data \$file ${summ_insertsize_hg_file}
        done
    echo -e "Sample\tinsert_size\tAll_Reads.fr_count" > ${summ_insertsize_dominant_virus_file}
    for file in ${insertsize_dominant_virus_files};do
        get_data \$file ${summ_insertsize_dominant_virus_file}
        done
    echo -e "Sample\tinsert_size\tAll_Reads.fr_count" > ${summ_insertsize_corrected_virus_file}
    for file in ${insertsize_corrected_virus_files};do
        get_data \$file ${summ_insertsize_corrected_virus_file}
        done
    """

    stub:
    """
    touch insertsize.hg.summ.txt
    touch insertsize.dominant.virus.summ.txt
    touch insertsize.corrected.virus.summ.txt
    """
}

