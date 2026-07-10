#!/bin/bash

mkdir -p /mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Fastq_Reverse/Sample_QC_reports

for file in /mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Fastq_Reverse/*.fastq.gz; do

	if [[ "$file" == *"Undetermined"* ]]; then
		continue
	fi

	fastqc -t 20 "$file" -o /mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Fastq_Reverse/Sample_QC_reports
done
