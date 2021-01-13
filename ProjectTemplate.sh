# Bioinformatics project template

# Instructions:
# CHANGE the PROJECT_DIR directory path in the script below:

PROJECT_DIR=/home/dcuser/dc_workshop/

# DATA:
DATA_DIR=${PROJECT_DIR}/data
DATA_RAW=${DATA_DIR}/untrimmed_fastq
DATA_TRIMMED=${DATA_DIR}/trimmed_fastq
DATA_SUBSET=${DATA_DIR}/trimmed_fastq_small
DATA_REFERENCE=${DATA_DIR}/ref_genome

# RESULTS:
RESULTS_DIR=${PROJECT_DIR}/results
RESULTS_FASTQ=${RESULTS_DIR}/fastqc_untrimmed_reads
RESULTS_SAM=${RESULTS_DIR}/sam
RESULTS_BAM=${RESULTS_DIR}/bam
RESULTS_BCF=${RESULTS_DIR}/bcf
RESULTS_VCF=${RESULTS_DIR}/vcf


# DOCS:
DOCS_DIR=${PROJECT_DIR}/docs

# SCRIPTS:
SCRIPTS_DIR=${PROJECT_DIR}/scripts

mkdir -p ${PROJECT_DIR} \
${DATA_DIR} \
${DATA_RAW} \
${DATA_TRIMMED} \
${DATA_SUBSET} \
${DATA_REFERENCE} \
${RESULTS_DIR} \
${RESULTS_FASTQ} \
${RESULTS_SAM} \
${RESULTS_BAM} \
${RESULTS_BCF} \
${RESULTS_VCF} \
${DOCS_DIR} \
${SCRIPTS_DIR} \