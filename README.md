# **D&Dseq: Docking and Deamination followed by sequencing** 

(https://www.landaulab.org)

![plot](./main/figs/DnD_workflow.png)

The symphony of gene expression consists of an orchestra of regulatory factors that are directed through multiple interconnected epigenetic signals, including chromatin accessibility and histone modifications, that affect transcription factor (TF) binding. These complex networks have been shown to be disrupted in cells during aging, disease or cancer. However, mechanistic profiling of these intricate pathways through single-cell analysis is mostly lacking, due to technical limitations of current methods for mapping DNA:protein interactions in single cells. This has resulted in a large gap in knowledge about where TFs or other chromatin remodelers bind to DNA, and how this binding is perturbed in pathological contexts in human tissues. To address this challenge, we have developed a versatile, high-throughput single-cell immuno-tethering DNA footprinting method that couples a species-specific antibody-binding nanobody to a cytosine base editing enzyme. Combined with single-cell ATAC-seq, this approach enables the profiling of even weak or transient factor binding to DNA in single-cells through identifying cytosine to uracil edits in target genomic regions following tagmentation and sequencing. Importantly, this Docking & Deamination followed by sequencing (D&D-seq) technique is easily incorporated into common single-cell multiomics workflows, allowing for multimodal analysis of gene regulation in single cells. We demonstrate the ability of D&D-seq to precisely profile CTCF and GATA1 binding in bulk as well as single-cell analyses of both cell lines and primary human cells. We further integrated D&D-seq with single-cell genotyping to assess the effect of IDH2 mutations on CTCF binding in human clonal hematopoiesis, identifying altered CTCF binding patterns in mutant cells. Together, the ability to directly measure TF or chromatin remodeler binding in vivo using high-throughput single-cell sequencing platforms will empower novel discoveries relating to chromatin and transcriptional regulation in human cells across physiological and disease contexts, at unprecedented scale and resolution.


## Computational pipeline for D&D analysis

### Extraction of D&D signal

To summarize genome edits introduced by DnD in scATAC-seq data, a bam file produced by CellRanger was split into each cell type based on barcode sequence and cell annotation using sinto. For bulk ATAC-seq data, bam files aligned to the genome using BWA-MEM2 were analyzed. After splitting bam files, the following six steps were performed to analyze DnD-mediated genomic variants. 

1. First, each bam file was preprocessed to remove uninformative and low quality read alignments. Duplicated reads were marked and simultaneously filtered using “picard MarkDuplicates” with a REMOVE_DUPLICATES=true parameter. Read alignments with high mapping quality Phred score (>= 20), primary alignment, reads aligned to intact chromosomes, and those with properly aligned mates were retained using samtools (version 1.19). 
2. Next, all single nucleotide variants (SNVs) found in each filtered bam file were collected using “bcftools mpileup” with following parameters, -a FORMAT/AD,FORMAT/DP,INFO/AD --no-BAQ --min-MQ 1 --max-depth 8000. The pileup result subsequently converted into the vcf format reporting SNVs supported by at least two mutant reads from minimum three aligned reads using an in-house Python script. 
3. Germline mutations were then filtered based on loci and alleles from the gnomAD database and variant allele frequency higher than 10%. If available, custom databases were provided to additionally filter uninformative mutations. 
4. Preprocessed bam files from step 1 were analyzed using MACS2 to call peaks with -f BAM --nomodel parameters. Peaks were then filtered with blacklist region annotation using bedtools (version 2.31.1). Motif analysis was performed using MEME Simple Enrichment Analysis (SEA) with HOmo sapiens COmprehensive MOdel COllection (HOCOMOCO) v11 core motif set to identify binding sites in peaks.
5. Target peaks were resized to 200 bp (up/downstream 100bp from the motif center) and overlaid with C-to-T and G-to-A variants identified in step 3. When multiple motif positions were found, a position with the highest score was chosen. Background peaks were resized to 200 bp by taking +/- 100bp from the peak summit. 


![plot](./main/figs/analytic_pipeline.png)

### Evaluation of D&D signal

To compare DnD edit counts between target and background peaks, normalized edit counts and signal-to-noise ratio (SNR) were calculated. 

First, SNVs counted in target and background peak groups were normalized by the total number of peaks and peak size and multiplied by a scaling factor of 100 to calculate the number of edits per 100 bp per peak. SNR was then calculated by dividing the DnD edit counts by the mean of non-DnD edit counts. 
Footprint analysis for DnD edits was performed by counting the number of DnD edits in each position from randomly sampled 200 target and background peaks. The random sampling was repeated for 10 times and mean and standard deviation were used for visualization.

![plot](./main/figs/DnD_signal.png)


### Prerequisite for D&D analysis

Following tools are required to run D&D analysis (more efficient version is coming). For most users, we recommend to install conda or mamba virtual environment to install these programs. For package installation and management, please advise with conda/mamba manuals.

| **tool** | **tested version** |
| -------- | -------- |
| Pytohn   | 3.3.9.19 |
| picard   | 3.1.1    |
| samtools | 1.19     |
| bcftools | 1.19     | 
| bedtools | 2.31.1   |
| bedops   | 2.4.41   | 
| vcf2bed  | 2.4.41   |
| macs2    | 2.2.9.1  |
| HOMER    | 4.11.1   |
| MEME     | 5.5.5    |

| **database** | **note** |
| -------- | ------- |
| ENCODE blacklist   | https://doi.org/10.1038/s41598-019-45839-z   |
| HOCOMOCO motif DB | included in MEME package    |
| gnomAD germline DB | https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38/    | 

**important note**
Bam files for samples, cell cluster or subclusters should be prepared in the following manner. For single-cell data, bam files can be separated by user-provided clusters using sinto. Please refer to the sinto's manual (https://timoast.github.io/sinto/).

```
$ tree bams

bams
├── ca46_ctcf
│   ├── ca46.5perc.bam
│   └── ca46.5perc.bam.bai
└── k562_gata1
    ├── k562.10perc.bam
    └── k562.10perc.bam.bai
```

### Running D&D analytic pipeline

D&D signals are collected and evaluated with three-step python scripts. By default, D&D edits will be called with motif analysis using MEME Simple Enrichment Analysis (SEA). Users can instead run HOMER2 with their own build or ChIP-seq as reference. 

**Step. 1: preprocessing, collecting and filtering variants, and first round peak calling**

```
$ ./dnd_pt1.py -h

usage:  [-h] -d DIR [-o OUTPUT] [--thread THREAD] [--start START] [--end END] [--mapq MAPQ] [--chrom] [--smt-other OTHER] [--se] [--count COUNT] [--alt ALT] [--fasta FASTA] [--gnomad GNOMAD] [--pass-gnomad]
        [--custom CUSTOM [CUSTOM ...]] [--vaf VAF] [--snv SNV] [--gsize GSIZE] [--opt OPT] [--blacklist BLACKLIST] [--pass-bklist]

optional arguments:
  -h, --help            show this help message and exit
  -d DIR, --Dir DIR     directory path
  -o OUTPUT, --Output OUTPUT
                        [Global] (*optional) output directory path
  --thread THREAD       [Global] (*optional) number of threads; default is 4
  --start START         [Global] (*optional) the first directory index
  --end END             [Global] (*optional) the last directory index
  --mapq MAPQ           [Step 1] (*optional) threshold for mapping quality; default is 20
  --chrom               [Step 1] (*optional) do not filter non-chromosome
  --smt-other OTHER     [Step 1] (*optional) other filter parameters for samtools in "~"
  --se                  [Step 1] (*optional) input bam is single-end
  --count COUNT         [Step 2] (*optional) minimum total counts; default is 3
  --alt ALT             [Step 2] (*optional) minimum ALT counts; default is 2
  --fasta FASTA         [Step 2] (*optional) genome fasta (indexed) used in alignment; e.g. cellranger/fasta/genome.fa
  --gnomad GNOMAD       [Step 3] (*optional) path to gnomAD vcf file
  --pass-gnomad         [Step 3] (*optional) do not run gnomAD filering
  --custom CUSTOM [CUSTOM ...]
                        [Step 3] (*optional) custom vcf file(s) to filter, multiple files are accepted
  --vaf VAF             [Step 3] (*optional) filter mutations frequent than this value; e.g. 10 for 10perc; default is 10
  --snv SNV             [Step 4] (*optoinal) SNV patterns to search; e.g. --snv "C>T,G>A, G>C"; default is "C>T,G>A"
  --gsize GSIZE         [Step 5] (*optional) effective genome size for macs2 callpeak; default is hs (homo sapiens)
  --opt OPT             [Step 5] (*optional) other parameters for macs2 callpeak
  --blacklist BLACKLIST
                        [Step 5] (*optional) blacklist file; default is hg38-blacklist.v2.bed
  --pass-bklist         [Step 5] (*optional) do not run blacklist filtering
```

```
./dnd_pt1.py -d <path_to_bam_directory> -o <path_to_output_directory> --thread 12 
```

Expected outputs are five directories in <-o>
```
├── step1_preprocess
│   ├── ca46_ctcf
│   │   ├── ca46.5perc.bam
│   │   ├── ca46.5perc.bam.bai
│   │   └── step1_picard_deduplication.txt
│   └── k562_gata1
│       ├── k562.10perc.bam
│       ├── k562.10perc.bam.bai
│       └── step1_picard_deduplication.txt
├── step2_mpileup
│   ├── ca46_ctcf
│   │   └── ca46.5perc.pileup
│   └── k562_gata1
│       └── k562.10perc.pileup
├── step3_fltvcf
│   ├── ca46_ctcf
│   │   └── ca46.5perc.flt.vcf
│   └── k562_gata1
│       └── k562.10perc.flt.vcf
├── step4_snvs
│   ├── ca46_ctcf
│   │   ├── ca46.5perc.flt.CT_GA.SNVs.vcf
│   │   ├── ca46.5perc.flt.SNVs.vcf
│   │   ├── ca46.5perc.snvs.bam
│   │   └── ca46.5perc.snvs.bam.bai
│   └── k562_gata1
│       ├── k562.10perc.flt.CT_GA.SNVs.vcf
│       ├── k562.10perc.flt.SNVs.vcf
│       ├── k562.10perc.snvs.bam
│       └── k562.10perc.snvs.bam.bai
└── step5_peaks
    ├── ca46_ctcf
    │   ├── peaks_bkflt.narrowPeak
    │   └── summits_bkflt.bed
    └── k562_gata1
        ├── peaks_bkflt.narrowPeak
        └── summits_bkflt.bed

```


**Step. 2: Joint peak calling and motif searching**

```
$ ./dnd_pt2.py -h

usage: [-h] -d DIR [-o OUTPUT] --mode [{sea,homer2}] [--gsize GSIZE] [--opt OPT] [--blacklist BLACKLIST] [--pass-bklist] [--motif MOTIF] [--homer-ref HM2REF]

optional arguments:
  -h, --help            show this help message and exit
  -d DIR, --Dir DIR     directory path
  -o OUTPUT, --Output OUTPUT
                        [Global] (*optional) output directory path
  --mode [{sea,homer2}]
                        [Global] which mode: "sea", "homer2"; default is "sea"
  --gsize GSIZE         [Step 5] (*optional) effective genome size for macs2 callpeak; default is hs
  --opt OPT             [Step 5] (*optional) other parameters for macs2 callpeak
  --blacklist BLACKLIST
                        [Step 5] (*optional) blacklist file; default is hg38-blacklist.v2.bed
  --pass-bklist         [Step 5] (*optional) do not run blacklist filtering
  --motif MOTIF         [Step 5, --mode:sea] (*optional) motif reference; default is <HOCOMOCOv11_core_HUMAN_mono_meme_format.meme>
  --homer-ref HM2REF    [Step 5, --mode:homer2] (*optional) homer2 reference"
```

```
./dnd_pt2.py -d <path_to_pt1_output_directory> --mode sea
```

This step will add a "<merged>" directory with joint peak calling results in <step5_peaks>, and intersected peaks in each sample's directory. Original MACS2 files will be stored in <celltype_specific> directory in each sample.

```
└── step5_peaks
    ├── ca46_ctcf
    │   ├── celltype_specific
    │   │   ├── peaks_bkflt.narrowPeak
    │   │   └── summits_bkflt.bed
    │   ├── peaks_bkflt.narrowPeak
    │   └── summits_bkflt.bed
    ├── k562_gata1
    │   ├── celltype_specific
    │   │   ├── peaks_bkflt.narrowPeak
    │   │   └── summits_bkflt.bed
    │   ├── peaks_bkflt.narrowPeak
    │   └── summits_bkflt.bed
    └── merged
        ├── meme_sea
        │   ├── sea.html
        │   ├── sea.tsv
        │   ├── sequences.tsv
        │   └── sites.tsv
        ├── peaks_bkflt.fasta
        ├── peaks_bkflt.narrowPeak
        ├── summits_bkflt.bed
        └── unfiltered
            ├── peaks_bkflt.narrowPeak
            └── summits_bkflt.bed
```

**Step. 3: Motif annotation and D&D edit evaluation**

```
$ ./dnd_pt3.py -h

usage: --mode only supports homer2 or sea [-h] -d DIR [-o OUTPUT] [--size SIZE] --sample SAMPLE [--var VARIANTS] --mode [{chip,homer2,sea}] [--chipseq CHIPSEQ] [--homer-ref HM2REF] [--motif MOTIF [MOTIF ...]]

optional arguments:
  -h, --help            show this help message and exit
  -d DIR, --Dir DIR     directory path
  -o OUTPUT, --Output OUTPUT
                        [Global] (*optional) output directory path
  --size SIZE           [Global] (*optional) test peak width size; default is 200
  --sample SAMPLE       [Global] "sample name" or "all" for all samples in <step5>
  --var VARIANTS        [Global] (*optional) expected D&D variants; default is "C>T,G>A"
  --mode [{chip,homer2,sea}]
                        [Global] which mode: "chip", "homer2" or "sea"
  --chipseq CHIPSEQ     [Step 6, --mode:chip] chip-seq reference
  --homer-ref HM2REF    [Step 6, --mode:homer2] (*optional) homer2 reference"
  --motif MOTIF [MOTIF ...]
                        [Step 6, --mode:homer2 or sea] <path to the homer2 motif file> for "homer2" or <TF name> for "sea"
```



Raw data is available at Gene Expression Omnibus ([GSEXXXXXX](https://www.landaulab.org))

