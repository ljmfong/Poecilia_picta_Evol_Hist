# _Poecilia picta_ Evolutionary History Pipeline
This directory contains scripts and pipelines used in the Evolutionary history of the Poecilia picta sex chromosome by: xxxx.
For further questions, please contact Lydia Fong: ljmfong@zoology.ubc.ca

The following pipeline is split into # parts and is labelled in the respective directories.

1. _P. picta_ SEX-DETector (inferring sex-linked sequences)
2. _P. picta_ PAML (phylogenetic analysis)
3. Characterization of P. picta Y Chromosome
4. Identiftying the ancestral sequence of P. picta
--------------------------------------------------------------------------------------------------------------------------------------------------------------------

## _P. picta_ SEX-DETector

Pipeline to analyse SNP segregation patterns and to identify sex-linked genes from RNA-seq data of a cross (parents and offspring of both sexes) using **[SEX-DETector](https://gitlab.in2p3.fr/sex-det-family/sex-detector/-/tree/master/)**

### a. Filter Fastq Files
**[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**: Quality check your fastq files. Example command:

* Important QC reports include: Per base sequence quality, Per sequence quality scores, Per seq GC content, Overrepresented sequences

**[MultiQC](https://multiqc.info/)**: Merge multiple FastQC files. Example command:

**[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)**: Trim off adaptors. Example command:

* Run FastQC on the trimmed sequences and compare results.

**[HISAT2](https://daehwankimlab.github.io/hisat2/)**: Map reads to verify sample integrity. Map trimmed reads to a reference genome (Poecilia reticulata) and check the mapping rate for each sample. 

* Build HISAT2 Index: Build index from reference genome
* Run HISAT2 Mapping:


### b. Constructring _de novo_ assembly

**[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)**: Build _de novo_ assembly. Trinity has a three step-process (Inchworm, Chrysalis, and Butterfly), some of the steps can take a long time. It is recommended to build the assembly in steps.

* Make a tab-delimited file that includes all trimmed forward and reverse sample file names (together with their paths). This will ensure that all samples are used to create the assembly. File should be of the format: 
sample_name \t sample_forward.fq.gz \t sample_reverse.fq.gz 

Then run:
*  Inchworm step:
*  Chrysalis step:
*  Butterfly step:

* Note the Trinity output format. You will get a '**[Trinity.fasta](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Output-of-Trinity-Assembly)**' output in a new 'trinity_out_dir/' that is created when Trinity is run. Trinity clusters transcripts loosely into 'genes' and contains 'isoforms' of those genes.

Custom script is used to determine assembly statistics (**03.assembly_stats.py**).
* Script takes the trinity output folder, searches for the assembly fasta file and calculates basic stats.



### c. Filtering Transcriptome 
Steps to filter transcriptome, removing redundancy,non-coding RNA and transcripts without an open reading frame

* **[RSEM](https://deweylab.github.io/RSEM/)** & **[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) Using Trinity Scripts**: Map reads to the _de novo_ transcriptome. 
Trinity has scripts that measure gene expression (RSEM) of reads and aligns them to the _de novo_ assembly (Bowtie2). These scripts can be found in the 'util' directory when Trinity is installed. However, both softwards must be installed and their PATH must be set in order to run the Trinity scripts. Run the following scripts:
* Prepare reference for alignment and abundance estimation:
* Run alignment and abundance esetimation:

Custom scripts are used to filter transcriptome to remove redundancy (**04.get-best-isoform.py** and **05.get-best-isoform-fasta.py**):

* **[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)**: Identify non-coding RNA (ncRNA): Download fasta file of ncRNA from closely related species from **[ENSEMBL](http://uswest.ensembl.org/index.html)** e.g. Oryzias_latipes.MEDAKA1.ncrna.fa. If using BLAST on command line, **[BLAST Commandline Manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/)**.

* Build BLAST index:
* BLAST transcript to ncRNA of related species:

Custom scripts are used to identify top BLAST hit (**06.get-ncrna.py**) and filter the assembly (**07.filter-assembly-ncrna.py**):

* **[Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki)**: Remove transcripts without open reading frames (ORFs).
* Extract long ORFs:
* Final coding region predictions:


### d. SEX-DETector
SEX-DETector requires a gametologs to coassemble into one single transcript. We need to further assemble the transcripts into contigs.

* **[CAP3](http://seq.cs.iastate.edu/cap3.html)**: Assemble contigs.
* Merge CAP3 singlets and contigs:

* **[BWA](http://bio-bwa.sourceforge.net/)**: Map trimmed reads to the final assembly
* Build BWA index:
* Align reads:
* Generate SAM format:

* **[SAMtools](http://www.htslib.org/doc/samtools.html)**: Convert SAMs to BAMs
* Build index:
* Compress into BAM files:
* Order individual BAm files:

* **[reads2snp](https://kimura.univ-montp2.fr/PopPhyl/index.php?section=tools)**: Genotype individuals at each locus. Note the output files:
* .alr - 
* .gen -
* Generate the gen_summary file -
* Run SEX-DETector


--------------------------------------------------------------------------------------------------------------------------------------------------------------------

## _P. picta_ PAML
