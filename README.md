# Poecilia_picta_Evol_Hist
This directory contains scripts and pipelines used in the Evolutionary history of the Poecilia picta sex chromosome by: xxxx.
For further questions, please contact Lydia Fong: ljmfong@zoology.ubc.ca

The following pipeline is split into # parts and is labelled in the respective directories.

1. Poecilia picta SEX-DETector (inferring sex-linked sequences)
2. P. picta PAML (phylogenetic analysis)
3. Characterization of P. picta Y Chromosome
4. Identiftying the ancestral sequence of P. picta
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

## P. picta SEX-DETector

Pipeline to analyse SNP segregation patterns and to identify sex-linked genes from RNA-seq data of a cross (parents and offspring of both sexes) using **[SEX-DETector](https://gitlab.in2p3.fr/sex-det-family/sex-detector/-/tree/master/)**

#### a. Filter Fastq Files
**[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**: Quality check your fastq files. Example command:

* Important QC reports include: Per base sequence quality, Per sequence quality scores, Per seq GC content, Overrepresented sequences

**[MultiQC](https://multiqc.info/)**: Merge multiple FastQC files. Example command:

**[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)**: Trim off adaptors. Example command:

* Run FastQC on the trimmed sequences and compare results.

**[HISAT2](https://daehwankimlab.github.io/hisat2/)**: Map reads to verify sample integrity. Map trimmed reads to a reference genome (Poecilia reticulata) and check the mapping rate for each sample. 

* Build HISAT2 Index: Build index from reference genome
* Run HISAT2 Mapping:


#### b. Constructring _de novo_ assembly

**[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)**: Build _de novo_ assembly. Trinity has a three step-process (Inchworm, Chrysalis, and Butterfly), some of the steps can take a long time. It is recommended to build the assembly in steps.

* Make a tab-delimited file that includes all trimmed forward and reverse sample file names (together with their paths). File should be of the format: 
sample_name \t sample_forward.fq.gz \t sample_reverse.fq.gz 

Then run:
*  Inchworm step:
*  Chrysalis step:
*  Butterfly step:



