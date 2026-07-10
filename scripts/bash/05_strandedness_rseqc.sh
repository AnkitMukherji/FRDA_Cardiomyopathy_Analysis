#!/bin/bash

INPUT="/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Dragen_Mapping_Fastq_Reverse_Trimmed"
REF="/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Strandedness/gencode.v49.annotation.bed12"
OUTPUT="/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Strandedness/rseqc_output"

for bam_file in "$INPUT"/*.bam;
do
	base_name=$(basename "$bam_file" .bam)
	infer_experiment.py -r "$REF" -i "$bam_file" > "${OUTPUT}/${base_name}_strandedness.txt"
done
