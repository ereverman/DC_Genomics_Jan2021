# Assessing Read Quality

## Questions:
* How can I describe the quality of my data?

## Objectives:
* Understand hoe a FASTQ file encodes per-base quaity scores.
* Interpret a FastQC plot summarizing per-base quality across all reads.
* Use for loops to automate operations on multiple files.

## Bioinformatic Workflows:
* The raw reads you get off the sequencer will need to be passed through different tools that assess and filter based on quality.
* Things that can affect read quality:
  * The length of the read (base quality decreases as a function of length)
  * Sequencing errors (dependent on DNA quality, sequencing technology)
  * Upstream sample processing (may use a custom kit, tissue used has unique difficulties)
* The execution of the tools to prepare raw sequencing files for downstream analyses is referred to as a pipeline or workflow
* Well-designed pipelines have the advantage of producing the correct input files for the next step in the pipeline
  * Files are generated according to a common standard format
* Pipelines can be fully automated once the bugs are worked out
* The backbone of the pipeline we will use can serve as a template or starting point for other pipelines you may need
* There are other tools available for checking quality and performing alignments
  * Keep in mind you should always assess the suitability and options of each tool for your own analyses
  
## Our Pipeline Overview:
* You will start with raw reads in FASTQ file format as if you received them from the sequencer
* We will perform quality assessments and filtering, which generates filtered FASTQ files
* We will then align the reads to the genome which will step through file conversions from FASTQ to SAM to BAM file formats
* We will perform post-alignment clean up with BAM files
* We will perform variant calling, which generates VCF files
  * At this point, we will be able to visualize variants and examine how they relate to our scientific question.
  
# Questions about the overview?

## Our data files:
* We will work with 3 samples from the LTEE E. coli experiment
  * Generations 5000, 15000, and 50000
  * Expect that the population has evolved over this time, and we will examine these genetic changes with variant calling
* Data are paired-end, so there will be two files for every sample
  * Enables both ends of the DNA to be sequenced
  * Helps align reads over repetitive regions of the genome more precisely
    1. Sequencing starts with library preparation. Genomic DNA (or cDNA if starting from RNA) is fragmented size-selected. 5’ and 3’ adapters are ligated. (Adapters can include sample-specific indices).
    2. Library is loaded into a flow cell and the fragments bind to short sequences that are complementary to the adapters. Bridge amplification generates clusters of primed templates across the flow cell.
    3. Sequencing by synthesis occurs as bases are added to the sequencing primer in a series of cycles. Bases are fluorescently labeled and when they are incorporated, they emit a base-specific wavelength that is decoded in an output file. One base is incorporated each cycle. The number of cycles determines the length of the read.
  * For a good resource: https://www.illumina.com/content/dam/illumina-marketing/documents/products/illumina_sequencing_introduction.pdf

* We are using the European Nucleotide Archive to access our data.
 * The ENA "provides a comprehensive record of the world's nucleotide sequencing information, covering raw sequencing data, sequence assemply information and functional annotation."
 * Data are provided in FASTQ format

## Download the data:
* All of the data files you've worked with so far were pre-loaded
* Under normal situations you will be up and downloading data
* These data may come from collaborators or from online repositories

* There are two programs that will download data from a remote server to your local or remote computer
 * wget = world wide web get, downloads web pages or data at a web address
 * curl = "see URL", displays webpages or data at a web address
* Typically use depends on operating system
* First determine whether you have wget or curl installed to use:
```
which curl
which wget

# Which is a BASH program that looks through everything you have installed and tells you where it is installed.
# no results means the searched program isn't installed
# both are installed on the AMI, try on your local computer
```

mkdir -p ~/dc_workshop/data/untrimmed_fastq/
cd ~/dc_workshop/data/untrimmed_fastq/

# -p option allows mkdir to create a new directory even if one of the parent direcoties doesn't already exist and supresses errors if the directory already exists WITHOUT overwriting that directory

curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz 

# This takes about 5 minutes for me.
# If there are connection issues, run this:
cd ~/dc_workshop/data/untrimmed_fastq/
cp ~/.backup/untrimmed_fastq/*fastq.gz .

```


## Unzip the data:
```
# The data files are very large and download in a compressed format.
ls -lh

# unzip:
gunzip SRR2584863_1.fastq.gz
ls -lh
# File is much bigger.
```

## Quality Control:
* FASTQ files have a specific format:
 * LINE 1 @ followed by information about the read
 * LINE 2 DNA sequence
 * LINE 3 + and sometimes same info as line 1
 * LINE 4 Quality scores for each base
 
```
head -n 4 SRR2584863_1.fastq
```
* Line 4 shows the quality for each nucleotide in the read.
* Quality is interpreted as the probability of an incorrect base call or base call accuracy.
* Numerical score is converted into a code so that the DNA sequence and the quality for each base match up perfectly.
 * Depends on the sequencing platform that generated the reads
 * Our data uses the standard Sanger quality PHRED score encoding, Illumina v 1.8
 * Quality scores range from 0 to 41 and relates to accuracy probability on a logarithmic scale.

```
# Looking back at the read, the last several bases have very poor quality scores (# = 2)
# Now, take a look at the last read in the file:

tail -n 4 SRR2584863_1.fastq

# Has a range of quality scores, seems to be better quality at the end of the read compared to the first read.
```
## Verify FASTQC is installed:
```
fastqc -h
# opens the manual page
# THUMBS UP/DOWN
```
## Running FASTQC:
```
cd untrimmed_fastq/

# There are a few ways to do this:
# Let's use some wildcards to simplify

fastqc *.fastq*

ls -lh
```
* FASTQC generates two types of files:
 * .zip = compressed set of multiple output files.
 * .html = stable webpage with a summary report for each sample
* Keep the raw data and results files separate:
```
mkdir -p ~/dc_workshop/results/fastqc_untrimmed_reads
mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads
mv *.htmp ~/dc_workshop/results/fastqc_untrimmed_reads
ls

cd ~/dc_workshop/results/fastqc_untrimmed_reads/
ls
```
## View FASTQC Results:
* To view the summary results, we need to access the html files
* Easiest for us to transfer them to our local computer from AWS using scp

```
# OPEN NEW TERMINAL TAB
mkdir -p ~/Desktop/fastqc_html

scp dcuser@XXX:~/dc_workshop/results/fastqc_untrimmed_reads/*.html ~/Desktop/fastqc_html
# Enter password

# ISSUES:
# macs, make sure you are running bash not zsh
cshs -s /bin/bash
# Enter password
```

## Do this:
* Open your file manager program and navigate to your fastqc directory on your desktop.
* Double click the first file.
* THUMBS UP/DOWN to verify everyone has access.

## FASTQC Results:
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/
* The first file:
 * A linked navigation summary is on the left
  * <ins>Basic Statistics</ins>
  * <ins>Per base sequence quality</ins>: Note that the axis isn't uniform; the later bins are aggregates of 5bp windows
  * <ins>Per tile sequence quality</ins>: Specific to Illumina library technology. Shows quality scores across all bases to see if there was a loss in quality associated with only one part of the flowcell (Could be caused by bubbles, smudges on flowcell, debrid in flowcell lane)
  * <ins>Per sequence quality scores</ins>: density plot of quality for all reads at all positions.
  * <ins>Per base sequence content</ins>: proportion of each base position over all of the reads. Expect ~25% through most of the read
  * <ins>Per sequence GC content</ins>: density plot of average GC content in each of the reads based on expected GC content. Deviation may be expected (especially for RNAseq data). Other issues may arise from contamination or library prep problems.
  * <ins>Per base N content</ins>: percent that 'N' occurs at a position in all reads. An increase at a particular position may indicate a problem occurred during sequencing.
  * <ins>Sequence length distribution</ins>: raw data, usually a sharp peak; trimmed data, may be broader distribution because reads are trimmed at different places.
  * <ins>Sequence duplication levels</ins>: expect most reads to occur only once. If sequences occur more than once, may indicate enrichment bias during PCR in library prep, or if sequencing was done at high coverage (eg RNA-seq, amplicon sequencing)
  * <ins>Overrepresented sequences</ins>: A list of sequences that occur more frequently than would be expected by chance.
  * <ins>Adapter content</ins>: when using long read lengths it is possible that some of the library inserts are shorter than the read length resulting in read-through to the adapter at the 3' end.
  * <ins>k-mer content</ins>: shows sequences which may show positional bias within the reads

## Breakout rooms:
* Breakout rooms for 5 minutes to look over each file and decide which samples look better/worse than others.


General patterns:
* Quality decreases toward the end of the read
* It is typical for read 1 scores to be higher than read 2
* Each of the samples has usable sequence, but the point varies among them for where they should be trimmed.


## FASTQC text output:
* Take a look at the .zip file.
* Go back to AWS and to the fastqc_untrimmed_reads/ directory
```
cd ~/dc_workshop/results/fastqc_untrimmed_reads/ 
ls 
```
* The .zip files are compressed files
* Contain multiple output files for each FASTQ file
* Use unzip to decompress the files:
```
unzip *.zip
ls -lh
unzip -h

# Wildcards don't work with unzip as expected because unzip expects only one file
```
* We could unzip each file one by one
* What happens when you have to do that for 500 files?
* Let's write a for loop instead.

```
for filename in *.zip
do
unzip $filename
done

# Using some symbols here that have different meanings in different contexts
# $ is typically a prompt, but here we are using it as a symbol to assign a variable
# > is used as a redirect to send file output somewhere we specify but here it is a modified prompt

ls -F
tree
# shows that each of the directories has several types of contents
```
* unzip creates a new directory with subdirectories for each of the 6 samples
* We will focus on the summary .txt file

```
ls -F SRR2584863_1_fastqc/ 
less SRR2584863_1_fastqc/summary.txt

# Looks familiar, right?
```
## Documenting our work:
* Individual summaries are useful, but better to have them all together for a project summary
```
mkdir -p ~/dc_workshop/docs/
cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt

cd ~/dc_workshop/docs/
less fastqc_summaries.txt
```

## Challenge: Which samples failed at least one of FASTQC's quality tests? What tests did those samples fail?
* Hint: think back to shell lessons.

```
grep FAIL fastqc_summaries.txt
```

## Building our pipeline:
