#!/bin/bash

BAM_FILES=/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Dragen_Mapping_Fastq_Reverse_Trimmed
DIR_REF=/mnt/faruq2/lab_data/reference/human/hg38/star_gtf/gencode49/gencode.v49.annotation.gtf
DIR_OUT=/mnt/faruq2/lab_users/ankit/rna_seq_asangla_cardio_12_06_26_lane_001/RNA_Feature_Counts

mkdir -p $DIR_OUT
mapfile -t BAMS < <(find "$BAM_FILES" -maxdepth 1 -name "*.bam" | sort)
echo "Found ${#BAMS[@]} BAM files."

THREADS=20               # Accelerate execution across 16 logical CPU cores (-T)
STRANDEDNESS=2           # Enforce reverse-stranded logic for FAST-1 deconvolution (-s); # 0 = unstranded, 1 = forward, 2 = reverse
FEATURE_TYPE="exon"      # Base mathematical mapping strictly on exonic overlaps (-t)
META_FEATURE="gene_id"   # Aggregate all exon counts to the overarching gene level (-g)

# The Core Summarization Algorithm Execution
#    -p                : Activates fragment-level counting for paired-end datasets.
#    --countReadPairs  : Count read pairs (fragments) instead of reads.
#    -M                : Forces the retention and analysis of multi-mapping fragments.
#    --fraction        : Prevents library inflation by distributing multi-mapped counts.
#    -B                : Require both ends mapped
#    -C                : Exclude chimeric pairs
#    -s 2              : Mandates directional strand verification.
#    -T 40             : Dictates the parallel processing bandwidth.
#    -t exon           : Restricts spatial boundary analysis to exons.
#    -g gene_id        : Establishes the final output matrix aggregation level.
#    -a                : Specifies the complex genomic topology reference file.
#    -o                : Defines the single unified output matrix destination.

featureCounts \
    -T "$THREADS" \
    -p \
    --countReadPairs \
    -B \
    -C \
    -s "$STRANDEDNESS" \
    -t "$FEATURE_TYPE" \
    -g "$META_FEATURE" \
    -a "$DIR_REF" \
    -o "$DIR_OUT/Cardiomyocyte_RNA_counts.txt" \
    "${BAMS[@]}"

echo "==========================================="
echo "featureCounts completed successfully."
echo "Output:"
echo "  Counts : $DIR_OUT/Cardiomyocyte_RNA_counts.txt"
echo "  Summary: $DIR_OUT/Cardiomyocyte_RNA_counts.txt.summary"
echo "==========================================="

# Utilize the 'tail' command to skip the first line containing the execution command header.
# Utilize the 'cut' command to extract Column 1 (GeneID) and all columns from 7 onwards (Sample Counts).
# Removing the directory paths and the .bam extension

tail -n +2 "$DIR_OUT/Cardiomyocyte_RNA_counts.txt" | cut -f 1,7- | awk 'BEGIN{FS=OFS="\t"} {for(i=1;i<=NF;i++) {sub(/^.*\//,"",$i); sub(/\.bam$/,"",$i)}} {print}' > "$DIR_OUT/Cardiomyocyte_RNA_counts_clean_matrix.txt"
