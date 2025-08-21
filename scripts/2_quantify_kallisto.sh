#!/bin/bash

# Kallisto quantification for paired-end reads

INDEX=~/Desktop/ABioinformatics/assignment/data/ref/mm10_transcripts.idx
BASE_DIR=~/Desktop/ABioinformatics/assignment/data
TRIM_DIR=$BASE_DIR/trimmed_reads
OUT_DIR=$BASE_DIR/kallisto_output

mkdir -p $OUT_DIR

SAMPLES=(SRR9685614 SRR9685616 SRR9685619 SRR9685620)

for SAMPLE in "${SAMPLES[@]}"; do
  echo "Quantifying $SAMPLE..."

  kallisto quant -i $INDEX \
    -o $OUT_DIR/$SAMPLE \
    -b 100 -t 4 \
    --plaintext \
    $TRIM_DIR/${SAMPLE}_trimmed_1.fastq.gz \
    $TRIM_DIR/${SAMPLE}_trimmed_2.fastq.gz
done

echo "✅ All samples quantified."

