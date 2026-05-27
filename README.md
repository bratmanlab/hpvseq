<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/bratmanlab-hpvseq_logo_dark.png">
    <img alt="bratmanlab/hpvseq" src="docs/images/bratmanlab-hpvseq_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/bratmanlab/hpvseq)
[![GitHub Actions CI Status](https://github.com/bratmanlab/hpvseq/actions/workflows/nf-test.yml/badge.svg)](https://github.com/bratmanlab/hpvseq/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/bratmanlab/hpvseq/actions/workflows/linting.yml/badge.svg)](https://github.com/bratmanlab/hpvseq/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/hpvseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![bratmanlab template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/bratmanlab/tools/releases/tag/3.5.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/bratmanlab/hpvseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23hpvseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/hpvseq)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/bratmanlab)

## Introduction

**bratmanlab/hpvseq** is a bioinformatics nextflow pipeline including HPV DNA detection, HPV genotyping and quantification, integration detection, fragment insert size estimation and variant calling on **dual-UMIs** and **paired-end barcoding sequencing** data by demultiplexing method of ConsensusCruncher (https://github.com/pughlab/ConsensusCruncher). It is tested on HPC, and will be on cloud in future. 

<!-- TODO bratmanlab:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to bratmanlab here, in 15-20 seconds. For an example, see https://github.com/bratmanlab/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO bratmanlab: Include a figure that guides the user through the major workflow steps. Many bratmanlab
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->
<!-- TODO bratmanlab: Fill in short bullet-pointed list of the default steps in the pipeline -->1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and bratmanlab, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

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
nextflow run bratmanlab/hpvseq \
   -profile <hpc/../docker/singularity/.../institute> \
   -c hpcgenomes.config \
   -params-file cohort.yml ## input, outdir, ...
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/hpvseq/usage) and the [parameter documentation](https://nf-co.re/hpvseq/parameters).

## Pipeline 

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/hpvseq/results) tab on the bratmanlab website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/hpvseq/output).

```text
Flow
.   
├── fastq or SRA? ✅  
├── Pre-QC 
│   └── Sample index error rate matched? ✅ 
├── QC
│   ├── FastqQC 
│   └── MultiQC
│       └── [ Report: QC Sequencing ] ✅
├── Post-QC 
│   ├── Trim read length? [ included but not tested yet ] ✅ 
│   └── Merge fastq files? ✅ 
├── Alignment to human genome
│   ├── [ Report: QC coverage ] ✅
│   ├── [ Report: on-target rate ] ✅
│   ├── Dedup
│   │   └── SMaSh
│   │       └── [ Report: Sample swap ]  ✅ 
│   └── GATK BQSR  
│       └── ConsensusCruncher 
│           ├── dcs_sc: dcs.sc
│           │   └── [ Report: Variant calls ] ❓
│           ├── dcs_sc: all.unique.dcs
│           │   └── [ Report: Variants VAF ] ❓
│           ├── [ Report: Saturation rate ]  ✅
│           └── [ Report: QC depth ] ✅
├── Genotyping by aligning hg-unmapped reads to 38 HPV genomes
│   └── ConsensusCruncher 
│       └── dcs_sc: all.unique.dcs
│           └── [ Report: detected genotypes (-f 2 -q 30) ] ✅
├── Quantification preparation on dominant genotype
│   └── Alignment of hg-unmapped reads to the given genotype 
│       └── ConsensusCruncher 
│           ├── dcs_sc: dcs.sc
│           │   └── [ Report: Variant calls ] ❓
│           ├── dcs_sc: all.unique.dcs
│           │   └── [ Report: Variants VAF ] ❓
│           ├── [ Report: Saturation rate ] ✅
│           └── [ Report: Qantification coverage (-f 2) ] ✅ 
│ 
├── Quantitfication preparation on baseline-corrected genotype
│   └──Alignment of hg-unmapped reads to the given genotype
│      └── ConsensusCruncher
│           ├── dcs_sc: dcs.sc
│           │   └── [ Report: Variant calls ] ❓
│           ├── dcs_sc: all.unique.dcs
│           │   └── [ Report: Variants VAF ] ❓
│           ├── [ Report: Saturation rate ] ✅
│           └── [ Report: Quantification coverage (-f 2) ] ✅
│
├── HPV integration by SearcHPV (UMI-trimmed reads)
│   ├── if genotype corrected based on baseline
│   │   └── [ Report: Breakpoints on corrected genotype ] ✅
│   └── if no genotype correction
│       └── [ Report: Breakpoints on dominant genotype ] ✅
├── Insert size by picard (properly-pairead reads; circle not corrected)
│   ├── if genotype corrected based on baseline
│   │   └── [ Report: Insert size metrics on corrected genotype (-f 2) ] ✅ 
│   └── if no genotype correction
│       └── [ Report: Insert size metrics on dominant genotype (-f 2) ] ✅
└── Genome-level distribution of depth of coverage for virus
    └── [ Report: depth of coverage/bp ] ❓ 

Output
.
├── multiqc
│   └── multiqc_data
├── tag_to_header
├── bwa
├── consensus
│   ├── hg19 
│   ├── HPV_38genomes
│   ├── corrected
│   └── dominant 
├── report
│   ├── coverage
│   ├── genotyping
│   ├── hsmetrics
│   ├── insertsize
│   ├── integration
│   └── smash
├── summary
└── logs 


Check completeness and failures on output files and logs
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
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
