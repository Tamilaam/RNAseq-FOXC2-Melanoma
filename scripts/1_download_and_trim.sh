#!/bin/bash

# Set base directory
BASE_DIR=~/Desktop/ABioinformatics/assignment/data
cd "$BASE_DIR"

# === 1. Download SRA files ===
SAMPLES=(SRR9685614 SRR9685616 SRR9685619 SRR9685620)

for SAMPLE in "${SAMPLES[@]}"; do
  prefetch "$SAMPLE"
done

# === 2. Convert to FASTQ ===
for SAMPLE in "${SAMPLES[@]}"; do
  fasterq-dump "$SAMPLE"
done

# === 3. Compress FASTQ files ===
for SAMPLE in "${SAMPLES[@]}"; do
  gzip ${SAMPLE}_1.fastq
  gzip ${SAMPLE}_2.fastq
done

# === 4. Trimming with Trimmomatic ===
TRIMMOMATIC_JAR=~/Desktop/ABioinformatics/assignment/tools/Trimmomatic-0.39/trimmomatic-0.39.jar
ADAPTERS=~/Desktop/ABioinformatics/assignment/tools/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa
TRIM_DIR=$BASE_DIR/trimmed_reads

mkdir -p "$TRIM_DIR"

for SAMPLE in "${SAMPLES[@]}"; do
  java -jar "$TRIMMOMATIC_JAR" PE -threads 8 \
    ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz \
    $TRIM_DIR/${SAMPLE}_trimmed_1.fastq.gz $TRIM_DIR/${SAMPLE}_unpaired_1.fastq.gz \
    $TRIM_DIR/${SAMPLE}_trimmed_2.fastq.gz $TRIM_DIR/${SAMPLE}_unpaired_2.fastq.gz \
    ILLUMINACLIP:$ADAPTERS:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70
done
