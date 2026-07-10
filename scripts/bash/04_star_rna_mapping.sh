#!/bin/bash

IN="/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Fastq_Reverse_Trimmed"
OUT="/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Star_Mapping"

mkdir -p "$OUT"

for R1 in "$IN"/*R1.fastq.gz
do
    R2="${R1/_R1/_R2}"

    # Extract sample name
    sample=$(basename "${R1%%_trimmed_L001_R*}")

    # Check if mate pair exists
    if [[ ! -f "$R2" ]]; then
        echo "ERROR: Missing R2 file for $R1"
        continue
    fi

    TMP_DIR="$OUT/tmp_${sample}"

    echo "========================================="
    echo "Processing sample: $sample"
    echo "R1: $R1"
    echo "R2: $R2"
    echo "========================================="

    STAR \
        --runThreadN 20 \
        --genomeDir /mnt/faruq2/lab_users/ankit/lcl_rnaseq/ref_gtf/genome_index/ \
        --readFilesCommand zcat \
        --readFilesIn "$R1" "$R2" \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --outTmpDir "$TMP_DIR" \
        --outFileNamePrefix "$OUT/${sample}_"

    # Remove temporary files after successful completion
    rm -rf "$TMP_DIR"

done
