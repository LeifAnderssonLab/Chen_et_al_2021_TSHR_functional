# Commands to calculate pi and Tajima's D from Pool-seq data using PoPoolation 1.2.2
# Angela P. Fuentes-Pardo, email: apfuentesp@gmail.com
# Uppsala University
# Date: 2021-04-23

#!/bin/bash

# Load libraries (in Uppmax computer cluster)
module load bioinfo-tools
module load samtools/1.10

# Set environment variables
WORK_DIR='path-to-working-directory'
CHR='chr-of-interest'
REF_FASTA='path-to-reference-fasta'
BAM_FILE='path-to-bam-file'
SAMPLE_ID='sample-name'
SUFFIX='output-name-suffix'
POPOOLATION_v1_PATH='path-to-local-copy-Dir-of-popoolation1-scripts'
MIN_COV='min-coverage'
MAX_COV='max-coverage'
WINDOW_SIZE='sliding-window-size'
STEP_SIZE='step-size'
POOL_SIZE='ploidy*sample-size'

# Go to working directory
cd $WORK_DIR

# Create a mpileup file
echo Creating mpileup ...
samtools mpileup -B -Q 0 -q 20 -r ${CHR} -f ${REF_FASTA} ${BAM_FILE} -o ${SAMPLE_ID}.${CHR}.$SUFFIX.pileup

# Filter out indels and +- 5bp around
echo Filtering indels ...
perl ${POPOOLATION_v1_PATH}/basic-pipeline/identify-genomic-indel-regions.pl --indel-window 5 --min-count 2 --input ${SAMPLE_ID}.${CHR}.$SUFFIX.pileup --output ${SAMPLE_ID}.${CHR}.$SUFFIX.indels.gtf
perl ${POPOOLATION_v1_PATH}/basic-pipeline/filter-pileup-by-gtf.pl --input ${SAMPLE_ID}.${CHR}.$SUFFIX.pileup --gtf ${SAMPLE_ID}.${CHR}.$SUFFIX.indels.gtf --output ${SAMPLE_ID}.${CHR}.$SUFFIX.idl.pileup && rm ${SAMPLE_ID}.${CHR}.$SUFFIX.pileup

# Subsample to uniform coverage
echo Subsampling to uniform coverage ...
perl ${POPOOLATION_v1_PATH}/basic-pipeline/subsample-pileup.pl --min-qual 20 --method withoutreplace --max-coverage ${MAX_COV} --fastq-type sanger --target-coverage ${MIN_COV} --input ${SAMPLE_ID}.${CHR}.$SUFFIX.idl.pileup --output ${SAMPLE_ID}.${CHR}.$SUFFIX.idl.scov.pileup

# Calculate nucleotide diversity (pi) and Tajima's D
echo Calculating pi ...
perl ${POPOOLATION_v1_PATH}/Variance-sliding.pl --fastq-type sanger --measure pi --input ${SAMPLE_ID}.${CHR}.$SUFFIX.idl.scov.pileup --min-count 2 --min-coverage ${MIN_COV} --max-coverage ${MAX_COV} --min-covered-fraction 0.4 --pool-size ${POOL_SIZE} --window-size ${WINDOW_SIZE} --step-size ${STEP_SIZE} --output ${SAMPLE_ID}.${CHR}.$SUFFIX.idl.scov.pi.txt #--snp-output ${SAMPLE_ID}.${CHR}.$SUFFIX.idl.scov.snps

echo Calculating Tajimas D ...
perl ${POPOOLATION_v1_PATH}/Variance-sliding.pl --fastq-type sanger --measure D --input ${SAMPLE_ID}.${CHR}.$SUFFIX.idl.scov.pileup --min-count 2 --min-coverage ${MIN_COV} --max-coverage ${MAX_COV} --min-covered-fraction 0.4 --pool-size ${POOL_SIZE} --window-size ${WINDOW_SIZE} --step-size ${STEP_SIZE} --output ${SAMPLE_ID}.${CHR}.$SUFFIX.idl.scov.TajD.txt
