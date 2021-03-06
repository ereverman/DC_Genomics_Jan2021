# Automating a Variant Calling Workflow

## Objectives:
* Make a pipeline more efficient and less error-prone
* Write a shell script with multiple variables
* Incorporate a for loop into a shell script

## What is a shell script?
* We've written short shell scripts already
* Shell scripts can be many lines long and can be written to accomplish several steps
* Especially useful for bioinformatics pipelines

* An important element of our shell script is the for loop
  * We've written a few for loops already that will be useful
  
## Analyzing Quality with FastQC:
* Use touch to create a new file in a new scripts/ directory where we will write the shell script

```
mkdir -p ~/dc_workshop/scripts
cd ~/dc_workshop/scripts
touch read_qc.sh

ls -lh

# now open the file:

nano read_qc.sh

# Type the following into your shell script:

set -e  # ensures script exits if error occurs
cd ~/dc_workshop/data/untrimmed_fastq/  # move to our raw data directory

echo "Running FastQC..."  # status message
fastqc *.fastq*   # run fastqc on all files with .fastq* extension

echo "Saving FastQC results..."   # status message
mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads/    # Move results to new EXISTING directory
mv *.html ~/dc_workshop/results/fastqc_untrimmed_reads/

cd ~/dc_workshop/results/fastqc_untrimmed_reads/    # move to fastqc results directory

echo "Unzipping fastqc results..."    # status message for for loop
for filename in *.zip
do
unzip $filename
done

echo "Saving summary..."
cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt
```

* Let's see how it does:
```
bash read_qc.sh

# reports status messages
# asks if we want to replace files (because we already ran them)
# select A

# Check on our files
ls -lh ../results/fastqc_untrimmed_reads/
ls -lh ../docs/

```

## Automating Variant Calling:
* Similar process, but more steps and more complicated loops.
* Instead of typing the script out together, we will download a script

```
cd ~/dc_workshop/scripts
curl -O https://raw.githubusercontent.com/datacarpentry/wrangling-genomics/gh-pages/files/run_variant_calling.sh
ls

less run_variant_calling.sh
```
* change directory
* set a variable "Genome" to describe where to find the reference
* index the genome
* create a subdirectory structure
* for loop that iterates through each of the sample files to perform variant calling

```
bash run_variant_calling.sh
```

## So can we answer the biological question?
* How do the three samples (generation 5000, 15000, 50000) change over time?

```
ls ../results/vcf/

for infile in ~/dc_workshop/results/vcf/*_final_variants.vcf
do
echo ${infile}
grep -v "#" ${infile} | wc -l
done

# Sample 9044 generation 5000 = 10 mutations
# Sample 4863 generation 15000 = 25 mutations
# Sample 4866 generation 50000 = 766 mutations

# It is likely that a hupermutable phenotype evolved between generation 15000 and 50000
```

## We are done, but here is a challenge homework:
* Our alignment was performed using a subset of the data.
* Copy and edit the variant calling script to run on the full data and see how the mutations differ

```
cd ../data/trimmed_fastq/
gunzip *.trim.fastq.gz

# edit variant calling script

set -e
cd ~/dc_workshop/results

genome=~/dc_workshop/data/ref_genome/ecoli_rel606.fasta

bwa index $genome

mkdir -p sam_full bam_full bcf_full vcf_full

for fq1 in ~/dc_workshop/data/trimmed_fastq/*_1.trim.fastq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.fastq)
    echo "base name is $base"

    fq1=~/dc_workshop/data/trimmed_fastq/${base}_1.trim.fastq
    fq2=~/dc_workshop/data/trimmed_fastq/${base}_2.trim.fastq
    sam=~/dc_workshop/results/sam_full/${base}.aligned.sam
    bam=~/dc_workshop/results/bam_full/${base}.aligned.bam
    sorted_bam=~/dc_workshop/results/bam_full/${base}.aligned.sorted.bam
    raw_bcf=~/dc_workshop/results/bcf_full/${base}_raw.bcf
    variants=~/dc_workshop/results/bcf_full/${base}_variants.vcf
    final_variants=~/dc_workshop/results/vcf_full/${base}_final_variants.vcf 

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam 
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
    done
    
    
for infile in ~/dc_workshop/results/vcf_full/*_final_variants.vcf
do
echo ${infile}
grep -v "#" ${infile} | wc -l
done

# Sample 9044 generation 5000 = 10 mutations
# Sample 4863 generation 15000 = 27 mutations
# Sample 4866 generation 50000 = 793 mutations
```

### Supporting reproducibility

```
touch sysinfo.txt
lshw > sysinfo.txt
apt list --installed –> pkginfo.txt

```
