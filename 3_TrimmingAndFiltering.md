# Trimming and Filtering

## Objectives:
* Remove sequence data that doesn't meet my quality standards with Trimmomatic.
* Select and set multiple options for command-line bioinformatic tools.
* Write for loops with two variables.


## Cleaning reads:
* FASTQC revealed that there is variability in sequencing quality
* We want to base our downstream analyses on data that has high quality and high confidence in accuracy
  * Low confidence base calls can inflate false positives and false negatives
  * Variant calling is sensitive to low confidence in base calls.

## Trimmomatic:
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
```
trimmomatic

THUMBS UP/DOWN
```
* We have to specify whether data is PE or SE
* There are a number of flags that allow you to optimize the tool for your system
* Next we have positional arguments -- order matters for how they are specified
* PE data: two input files, four output files
* Next, trimmomatic expects trimming parameters

```
cd ~/dc_workshop/data/untrimmed_fastq
ls -lh

# We have adapter sequences in our data files (saw with fastqc)
# Adapter sequences came with the installation of trimmomatic
# copy the adapter sequences to the current directory

cp ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/NesteraPE-PR.fa .

trimmomatic PE \
SRR2589044_1.fastq.gz SRR2589044_2.fastq.gz \
SRR2589044_1.trim.fastq.gz SRR2589044_1un.trim.fastq.gz \
SRR2589044_2.trim.fastq.gz SRR2589044_2un.trim.fastq.gz \
SLIDINGWINDOW:4:20 \		
MINLEN:25 \			
ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

```

* What percent of reads were discarded? 0.23%
* What percent of reads were kept in both pairs? 79.96%

* Confirm that we have the output files.
```
ls SRR25489044*

# trimmed files should be smaller (we've removed reads)
ls  SRR25489044* -lh
```
* We've trimmed our fastq files successfully!
* It was kind of a pain to type that out once.

```
# We want to loop through all the files. This is easiest if they are all the same type of file
gzip SRR2584863_1.fastq 

for infile in *_1.fastq.gz
do
  base=$(basename ${infile} _1.fastq.gz)
  trimmomatic PE ${infile} ${base}_2.fastq.gz \
               ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
               ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
               SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
done

```
