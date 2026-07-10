#!/bin/bash

GTF="/mnt/faruq2/lab_data/reference/human/hg38/star_gtf/gencode49/gencode.v49.annotation.gtf"
INPUT="/mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Dragen_Mapping_Fastq_Reverse_Trimmed"
OUTPUT_DIR="/mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Qualimap_Results"

mkdir -p "$OUTPUT_DIR"

for bam in "$INPUT"/*.bam
do
    sample=$(basename "$bam" .bam)

    sample_outdir="${OUTPUT_DIR}/${sample}"
    mkdir -p "$sample_outdir"

    export _JAVA_OPTIONS="-Djava.io.tmpdir=$sample_outdir"

    qualimap rnaseq \
        --java-mem-size=50G \
        -bam "$bam" \
        -gtf "$GTF" \
        -outdir "$sample_outdir" \
        -outfile "${sample}_qualimap" \
        -p strand-specific-reverse \
        -pe
done
