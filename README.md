#**D&Dseq: Docking and Deamination followed by sequencing** 

(https://www.landaulab.org)

![plot](./main/figs/DnD_workflow.png)

Abstract: blah blah blah~


Computational pipeline for D&D analysis

**Extraction of D&D signal**

To summarize genome edits introduced by DnD in scATAC-seq data, a bam file produced by CellRanger was split into each cell type based on barcode sequence and cell annotation using sinto. For bulk ATAC-seq data, bam files aligned to the genome using BWA-MEM2 were analyzed. After splitting bam files, the following six steps were performed to analyze DnD-mediated genomic variants. 

1. First, each bam file was preprocessed to remove uninformative and low quality read alignments. Duplicated reads were marked and simultaneously filtered using “picard MarkDuplicates” with a REMOVE_DUPLICATES=true parameter. Read alignments with high mapping quality Phred score (>= 20), primary alignment, reads aligned to intact chromosomes, and those with properly aligned mates were retained using samtools (version 1.19). 
2. Next, all single nucleotide variants (SNVs) found in each filtered bam file were collected using “bcftools mpileup” with following parameters, -a FORMAT/AD,FORMAT/DP,INFO/AD --no-BAQ --min-MQ 1 --max-depth 8000. The pileup result subsequently converted into the vcf format reporting SNVs supported by at least two mutant reads from minimum three aligned reads using an in-house Python script. 
3. Germline mutations were then filtered based on loci and alleles from the gnomAD database and variant allele frequency higher than 10%. If available, custom databases were provided to additionally filter uninformative mutations. 
4. Preprocessed bam files from step 1 were analyzed using MACS2 to call peaks with -f BAM --nomodel parameters. Peaks were then filtered with blacklist region annotation using bedtools (version 2.31.1). Motif analysis was performed using MEME Simple Enrichment Analysis (SEA) with HOmo sapiens COmprehensive MOdel COllection (HOCOMOCO) v11 core motif set to identify binding sites in peaks.
5. Target peaks were resized to 200 bp (up/downstream 100bp from the motif center) and overlaid with C-to-T and G-to-A variants identified in step 3. When multiple motif positions were found, a position with the highest score was chosen. Background peaks were resized to 200 bp by taking +/- 100bp from the peak summit. 


![plot](./main/figs/analytic_pipeline.png)

**Evaluation of D&D signal**

To compare DnD edit counts between target and background peaks, normalized edit counts and signal-to-noise ratio (SNR) were calculated. 

First, SNVs counted in target and background peak groups were normalized by the total number of peaks and peak size and multiplied by a scaling factor of 100 to calculate the number of edits per 100 bp per peak. SNR was then calculated by dividing the DnD edit counts by the mean of non-DnD edit counts. 
Footprint analysis for DnD edits was performed by counting the number of DnD edits in each position from randomly sampled 200 target and background peaks. The random sampling was repeated for 10 times and mean and standard deviation were used for visualization.

![plot](./main/figs/DnD_signal.png)

Raw data is available at Gene Expression Omnibus ([GSEXXXXXX](https://www.landaulab.org))

