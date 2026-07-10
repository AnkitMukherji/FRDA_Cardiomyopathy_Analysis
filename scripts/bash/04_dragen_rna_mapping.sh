#!/bin/bash

IN="/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Fastq_Reverse_Trimmed"
OUT="/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Dragen_Mapping_Fastq_Reverse_Trimmed"
REF="/mnt/faruq2/lab_data/reference/human/hg38/dragen_gtf"
GTF="/mnt/faruq2/lab_data/reference/human/hg38/dragen_gtf/gencode.v49.annotation.gtf"

for R1 in "$IN"/*R1.fastq.gz
do
    R2="${R1/_R1/_R2}"
    sample=$(basename "${R1%%_trimmed_L001_R*}")

    echo "Processing sample: $sample"
    echo "R1: $R1"
    echo "R2: $R2"

    dragen -f \
        -r "$REF" \
        -1 "$R1" \
        -2 "$R2" \
        -a "$GTF" \
        --RGID "$sample" \
        --RGSM "$sample" \
        --output-directory "$OUT" \
        --output-file-prefix "$sample" \
        --enable-rna true
done
