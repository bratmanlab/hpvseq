<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/bratmanlab-hpvseq_logo_dark.png">
    <img alt="bratmanlab/hpvseq" src="docs/images/bratmanlab-hpvseq_logo_light.png">
  </picture>
</h1>

**bratmanlab/hpvseq** is a bioinformatics nextflow pipeline including HPV DNA detection, HPV genotyping and quantification, and optional analyses including integration detection, fragment insert size estimation and variant calling on **dual-UMIs** and **paired-end barcoding sequencing** data by demultiplexing method of ConsensusCruncher (https://github.com/pughlab/ConsensusCruncher). It is tested on HPC, and will be on cloud in future. 

<!-- TODO bratmanlab:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to bratmanlab here, in 15-20 seconds. For an example, see https://github.com/bratmanlab/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO bratmanlab: Include a figure that guides the user through the major workflow steps. Many bratmanlab
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->

## Usage

<!-- TODO bratmanlab: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO bratmanlab: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run bratmanlab/hpvseq/main.nf \
   -profile <hpc/../docker/singularity/.../institute> \
   -c hpcgenomes.config \
   -params-file cohort.yml ## input, outdir, ...
```

## Pipeline 

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/hpvseq/results) tab on the bratmanlab website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/hpvseq/output).

> **Flow**
```text
.   
├── fastq or SRA?                                               [ SRA: storage+1 ] ✅  
├── Pre-QC 
│   └── Sample index error rate matched?                                           ✅ 
├─?─ QC
│   ├── FastqQC 
│   └── MultiQC
│       └── [ Report: QC Sequencing ]                                              ✅
├── Post-QC 
│   ├─?─ Trim read length?                                           [ storage+1 ] ✅ 
│   └── Merge fastq files?                                           [ storage+1 ] ✅ 
├── TagToHeader                                                      [ storage+1 ] ✅ 
├── Alignment to human genome                                        [ storage+1 ]
│   ├── [ Report: QC coverage ]                                                    ✅
│   ├── [ Report: on-target rate ]                                                 ✅
│   ├─?─ Dedup                                                      
│   │   └── SMaSh
│   │       └── [ Report: Sample swap ]                                            ✅ 
│   └── GATK BQSR                                                    [ storage+1 ] 
│       └── ConsensusCruncher                                        [ storage+1 ] 
│           ├── dcs_sc: dcs.sc
│           │   └── [ Report: Variant calls ]                                      ❓
│           ├── dcs_sc: all.unique.dcs
│           │   └── [ Report: Variants VAF ]                                       ❓
│           ├── [ Report: Saturation rate ]                                        ✅
│           └── [ Report: QC depth ]                                               ✅
├─?─ Genotyping by aligning hg-unmapped reads to 38 HPV genomes
│   └── ConsensusCruncher 
│       └── dcs_sc: all.unique.dcs
│           └── [ Report: detected genotypes (-f 2 -q 30) ]                        ✅
├── Quantification preparation on dominant genotype
│   └── Alignment of hg-unmapped reads to the given genotype 
│       └── ConsensusCruncher 
│           ├── dcs_sc: dcs.sc
│           │   └── [ Report: Variant calls ]                                      ❓
│           ├── dcs_sc: all.unique.dcs
│           │   └── [ Report: Variants VAF ]                                       ❓
│           ├── [ Report: Saturation rate ]                                        ✅
│           └── [ Report: Qantification coverage (-f 2) ]                          ✅ 
│ 
├─?─ Quantitfication preparation on baseline-corrected genotype
│   └──Alignment of hg-unmapped reads to the given genotype
│      └── ConsensusCruncher
│           ├── dcs_sc: dcs.sc
│           │   └── [ Report: Variant calls ]                                      ❓
│           ├── dcs_sc: all.unique.dcs
│           │   └── [ Report: Variants VAF ]                                       ❓
│           ├── [ Report: Saturation rate ]                                        ✅
│           └── [ Report: Quantification coverage (-f 2) ]                         ✅
│
├─?─ HPV integration by SearcHPV (UMI-trimmed reads)                  [ storage+1 ]
│   ├── if genotype corrected based on baseline
│   │   └── [ Report: Breakpoints on corrected genotype ]                          ✅
│   └── if no genotype correction
│       └── [ Report: Breakpoints on dominant genotype ]                           ✅
├─?─ Insert size (properly-pairead reads; circle not corrected)       
│   ├── if genotype corrected based on baseline
│   │   └── [ Report: Insert size metrics on corrected genotype (-f 2) ]           ✅ 
│   └── if no genotype correction
│       └── [ Report: Insert size metrics on dominant genotype (-f 2) ]            ✅
└─?─ Genome-level distribution of depth of coverage for virus
    └── [ Report: depth of coverage/bp ]                                           ❓
 
-?- represents that the process is skippable
* storage is roughly estimated on fastq files which accomodate temporary bam files 
```

> **Output**
```text
.
├── multiqc
│   └── multiqc_data
├── tag_to_header                                                    [ storage+1 ]
├── bwa                                                              [ storage+1 ]
├── consensus                                                        [ storage+1 ]
│   ├── hg19 
│   ├── HPV_38genomes
│   ├── corrected
│   └── dominant 
├── report
│   ├── coverage
│   ├── genotyping
│   ├── hsmetrics
│   ├── insertsize
│   ├── integration                                                  [ storage+1 ]
│   └── smash
├── summary
└── logs 

```

> **Check completeness and failures on output files and logs**
```text
.
├── Alignment (as basis)
│   └── bwa 
├── ConsensusCruncher ❓
│   ├── hg ❓
│   ├── genotyping ❓
│   ├── dominant genotypes ❓
│   └── corrected genotypes ❓
├── Integration ❓
│   └── on dominant or corrected genotypes ❓
├── Insert size ❓
│   ├── dominant genotypes ❓
│   └── corrected genotypes ❓

```

## Dependencies

SMaSH: https://github.com/rbundschuh/SMaSH

ConsensusCruncher: https://github.com/pughlab/ConsensusCruncher

SearcHPV: https://github.com/WenjinGudaisy/SearcHPV

## Credits

bratmanlab/hpvseq was originally written by Jinfeng Zou.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO bratmanlab: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#hpvseq` channel](https://nfcore.slack.com/channels/hpvseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO bratmanlab: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use bratmanlab/hpvseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO bratmanlab: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `bratmanlab` publication as follows:

> **The bratmanlab framework for community-curated bioinformatics pipelines.**
>
> Leung E, Han K, Zou J, Zhao Z, Zheng Y, Wang TT, Rostami A, Siu LL, Pugh TJ, Bratman SV. HPV Sequencing Facilitates Ultrasensitive Detection of HPV Circulating Tumor DNA. Clin Cancer Res. 2021 Nov 1;27(21):5857-5868. doi: 10.1158/1078-0432.CCR-19-2384. Epub 2021 Sep 27. PMID: 34580115; PMCID: PMC9401563.
> Han K, Zou J, Zhao Z, Baskurt Z, Zheng Y, Barnes E, Croke J, Ferguson SE, Fyles A, Gien L, Gladwish A, Lecavalier-Barsoum M, Lheureux S, Lukovic J, Mackay H, Marchand EL, Metser U, Milosevic M, Taggar AS, Bratman SV, Leung E. Clinical Validation of Human Papilloma Virus Circulating Tumor DNA for Early Detection of Residual Disease After Chemoradiation in Cervical Cancer. J Clin Oncol. 2024 Feb 1;42(4):431-440. doi: 10.1200/JCO.23.00954. Epub 2023 Nov 16. PMID: 37972346; PMCID: PMC10824379.
>
