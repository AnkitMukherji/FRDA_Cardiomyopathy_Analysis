#!/bin/bash

# -l denotes --length_required reads shorter than length_required will be discarded, default is 15. (int [=15])
# -p denotes --overrepresentation_analysis enable overrepresented sequence analysis.
# -5, --cut_front                      move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
# -3, --cut_tail                       move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
# -W, --cut_window_size                the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
# -M, --cut_mean_quality               the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])

mkdir -p \
/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Fastq_Reverse/trimmed_outputs \
/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Fastq_Reverse/qc_trimmed

OUT_TRIM=/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Fastq_Reverse/trimmed_outputs
OUT_QC_TRIM=/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Fastq_Reverse/qc_trimmed

for r1 in /mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Fastq_Reverse/*R1*.fastq.gz; do

    [[ "$r1" == *Undetermined* ]] && continue

    r2="${r1/_R1_/_R2_}"

    if [[ ! -f "$r2" ]]; then
        echo "Missing pair for $r1"
        continue
    fi

    sample=$(basename "${r1%%_L001_R*}")

    echo "Processing $sample"

    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$OUT_TRIM/${sample}_trimmed_L001_R1.fastq.gz" \
        -O "$OUT_TRIM/${sample}_trimmed_L001_R2.fastq.gz" \
        --detect_adapter_for_pe \
        -l 50 \
        -p \
        -5 -3 -W 4 -M 20 \
        --thread 20 \
        -h "$OUT_TRIM/${sample}_report.html" \
        -j "$OUT_TRIM/${sample}_report.json"

    echo "Running Trimmed QC..."

    fastqc \
        -t 20 \
        -o "$OUT_QC_TRIM" \
        "$OUT_TRIM/${sample}_trimmed_L001_R1.fastq.gz" \
        "$OUT_TRIM/${sample}_trimmed_L001_R2.fastq.gz"

done

echo "Workflow complete!"
