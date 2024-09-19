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


### Running D&D analytic pipeline

Following tools are required to run D&D analysis (more efficient version is coming). For most users, we recommend to install conda or mamba virtual environment to install these programs. For package installation and management, please advise with conda/mamba manuals.

  **tool**    **tested version**
  
  picard    3.1.1
  
  samtools    1.19
  
  bcftools    1.19
  
  bedtools    2.31.1
  
  bedops    2.4.41
  
  vcf2bed    2.4.41
  
  macs2    2.2.9.1
  
  HOMER    4.11.1
  
  MEME    5.5.5



Raw data is available at Gene Expression Omnibus ([GSEXXXXXX](https://www.landaulab.org))

