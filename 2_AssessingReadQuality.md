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

  
  * Paired-end reads can help with sequence assembly, especially if the reference genome is high quality (is with model organisms)


