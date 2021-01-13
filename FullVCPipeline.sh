# Full Variant Calling Sample Pipeline:

# INSTRUCTIONS:
# Edit path information to use with a new project (lines ).



set -e  # ensures script exits if error occurs



# STEP 1: FASTQC

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



# STEP 2: TRIM and FILTER

echo "Running Trimmomatic..."

cd ~/dc_workshop/data/untrimmed_fastq
cp ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/NexteraPE-PE.fa .

for infile in *_1.fastq.gz
do
base=$(basename ${infile} _1.fastq.gz)
trimmomatic PE ${infile} ${base}_2.fastq.gz \
${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
done

mv *.trim* ~/dc_workshop/data/trimmed_fastq

# STEP 3: ALIGN and PERFORM variant calling:


cd ~/dc_workshop/results

genome=~/dc_workshop/data/ref_genome/ecoli_rel606.fasta

bwa index $genome


for fq1 in ~/dc_workshop/data/trimmed_fastq_small/sub/*_1.trim.sub.fastq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.sub.fastq)
    echo "base name is $base"

    fq1=~/dc_workshop/data/trimmed_fastq_small/sub/${base}_1.trim.sub.fastq
    fq2=~/dc_workshop/data/trimmed_fastq_small/sub/${base}_2.trim.sub.fastq
    sam=~/dc_workshop/results/sam/${base}.aligned.sam
    bam=~/dc_workshop/results/bam/${base}.aligned.bam
    sorted_bam=~/dc_workshop/results/bam/${base}.aligned.sorted.bam
    raw_bcf=~/dc_workshop/results/bcf/${base}_raw.bcf
    variants=~/dc_workshop/results/bcf/${base}_variants.vcf
    final_variants=~/dc_workshop/results/vcf/${base}_final_variants.vcf 

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
    done