# Variant Calling Workflow

## Objectives
* Find sequence variants between my sample and a reference genome
* Describe types of data formats encountered during variant calling
* Use command line tools to perform variant calling

## Burrows Wheeler Aligner (BWA)
http://bio-bwa.sourceforge.net

* Two step process
  * Index reference genome
  * Align reads
  
## Setting up:
* First, need the reference genome
* Reference genomes are typically archived and downloadable
* Location varies by organism and by who generated and maintains the data
* References are also updated periodically

```
cd ~/dc_workshop
mkdir -p data/ref_genome
curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz

gunzip data/ref_genome/ecoli_rel606.fasta.gz

# Take a peek at the file:
head data/ref_genome/ecoli_rel606.fasta
# Chromosome name followed by the name of the sequence
```

* Variant calling can be a time consuming process
* To help us run the next commands, we will work with a subset of the data we generated

```
cd ~/dc_workshop

curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
ls -lh
tar xvf sub.tar.gz
ls -lh
mv sub/ ~/dc_workshop/data/trimmed_fastq_small

# Make directories to hold results files:
mkdir -p results/sam results/bam results/bcf results/vcf
```

## Index the reference genome:
```
bwa index data/ref_genome/ecoli_rel606.fasta
ls -lh data/ref_genome/ecoli_rel606.fasta

# Generated several files which will be used during the alignment
```

## Align reads to the reference:
* Using BWA-MEM which is optimized for longer reads generated on NGS technology
```
bwa mem data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam
```

## SAM/BAM format:
https://genome.sph.umich.edu/wiki/SAM
* tab-delimited txt file with information for each individual read and its alignment to the genome
* Sometimes contains a header with information about the source of data, alignment, reference sequence, etc.
* Contains 11 mandatory fields with information for mapping reads
* Additional fields may contain aligner-specific information
```
head results/sam/SRR2584866.aligned.sam
```
* BAM file: compressed binary version of the SAM file
  * Reduces size
  * Allows for indexing for efficient access of data

```
samtools view -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam

# -S indicates the input file is sam format and to output in -b bam format


ls -lh results/bam/
ls -lh results/sam/
```
## Sort BAM file by coordinates:
* The way that files are sorted differs between alignment tools
* Downstream tools require different types of sorted files
```
samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam

# -o indicates where to output the file

# Learn more about our new sorted bam file:
samtools flagstat results/bam/SRR2584866.aligned.sorted.bam
```

## Variant Calling:
* A variant is a different nucleotide at a given position compared to the reference
* Typically accompanied by an estimate of variant frequency and confidence
* Many tools available for variant calling
* We will use bcftools
https://samtools.github.io/bcftools/bcftools.html
* There are three steps

## Calculate read coverage:
* Describes the average number of reads that align to reference bases
* Determines whether variant discovery can be made with a certain degree of confidence
* Influenced by the size of the genome, read length, and the number of samples sequenced
* "Requirements" for coverage are influenced by experimental design and data type
```
bcftools mpileup -O b -o results/bcf/SRR2584866_raw.bcf \
-f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam 

# -O b tells bcftools to generate a bcf format output file
# -o specifies where to write the output
# -f specifies path to reference genome

# BONUS: What is the average coverage?
samtools view -H results/bam/SRR2584866.aligned.sorted.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'
# 4629812

samtools depth  results/bam/SRR2584866.aligned.sorted.bam  |  awk '{sum+=$3} END { print "Average = ",sum/4629812}'
# Average =  9.97389
```
## Detect SNPs
```
bcftools call --ploidy 1 -m -v -o results/bcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf

# --ploidy 1 for haploid organisms
# -m allows for multiallelic and rare-variant calling
# -v output variant sites only
# -o specifies where to write the output
```

## Filter and report SNPs in Variant Calling Format (VCF)
https://gatk.broadinstitute.org/hc/en-us
```
vcfutils.pl varFilter results/bcf/SRR2584866_variants.vcf  > results/vcf/SRR2584866_final_variants.vcf
less -S results/vcf/SRR2584866_final_variants.vcf
```

## Exercise:
How many variants are in the vcf file?

```
grep -v "#" results/vcf/SRR2584866_final_variants.vcf | wc -l

# 766
```

## Visualizing the alignment (two methods)
* How many have IGV installed on your local computer? THUMBS UP/DOWN

* Data visualization can give you ideas for further analyses and can help identify abnormalities
* We will look at two tools:

### Method one:
```
# First index the BAM file using samtools:
samtools index results/bam/SRR2584866.aligned.sorted.bam

# Then visualize mapped reads with tview:
samtools tview results/bam/SRR2584866.aligned.sorted.bam data/ref_genome/ecoli_rel606.fasta

# First line shows genome coordinates in reference
# Second line shows reference sequence
# Third line shows the concensus sequence determined from the sequence reads . is match to reference.
# Below the horizontal line we see reads aligned with the reference

# We can check a specific location:
g
CP000819.1:4377265
return

# G is the variant A is the reference
q to quit

```

### Method two:
http://software.broadinstitute.org/software/igv/download

* IGV is a GUI browser that can be installed on your local computer
* We need to transfer files to our local desktop

```
# Move to your local computer terminal
mkdir ~/Desktop/files_for_igv
cd ~/Desktop/files_for_igv

# We need to transfer 4 files:
scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/results/bam/SRR2584866.aligned.sorted.bam ~/Desktop/files_for_igv

scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/results/bam/SRR2584866.aligned.sorted.bam.bai ~/Desktop/files_for_igv

scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/data/ref_genome/ecoli_rel606.fasta ~/Desktop/files_for_igv

scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/results/vcf/SRR2584866_final_variants.vcf ~/Desktop/files_for_igv
```
1. Open IGV
2. Load reference genome file
 * Genomes --> Load Genomes from File
3. Load BAM file
 * File --> Load from File
4. Load VCF
 * File --> Load from File
