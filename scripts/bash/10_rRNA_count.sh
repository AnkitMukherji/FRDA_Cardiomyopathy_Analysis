#!/bin/bash

set -euo pipefail

###############################
# User settings
###############################

###############################
# User settings
###############################

IN="/mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Dragen_Mapping_Fastq_Reverse_Trimmed"
OUT="/mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Feature_Counts"
GTF_DIR="/mnt/faruq2/lab_data/reference/gencode/gencode49"
GTF="$GTF_DIR/gencode.v49.annotation.gtf"

RRNA_GTF="$GTF_DIR/gencode.v49.rRNA_only.gtf"

THREADS=20
STRAND=2

mkdir -p "$OUT"

LOG="$OUT/featureCounts_pipeline.log"

exec > >(tee -a "$LOG") 2>&1

echo "==========================================="
echo "RNA FeatureCounts + rRNA QC Pipeline"
date
echo "==========================================="

###############################
# Check input
###############################

echo
echo "Checking BAM files..."

BAMS=("$IN"/*.bam)

if [ ${#BAMS[@]} -eq 0 ]; then
    echo "ERROR: No BAM files found."
    exit 1
fi

echo "Found ${#BAMS[@]} BAM files."

###############################
# Create rRNA GTF
###############################

if [ ! -f "$RRNA_GTF" ]; then

    echo
    echo "Creating rRNA annotation..."

    awk '
    BEGIN{FS=OFS="\t"}

    /^#/{
        print
        next
    }

    $3=="exon" &&
    ($0~/gene_type "rRNA"/ || $0~/gene_type "Mt_rRNA"/)
    {
        print
    }
    ' "$GTF" > "$RRNA_GTF"

fi

###############################
# rRNA genes count
###############################

echo
echo "Checking rRNA annotation..."

NGENES=$(grep -v "^#" "$RRNA_GTF" | \
grep 'gene_id "' | \
sed 's/.*gene_id "\(.*\)".*/\1/' | \
sort -u | \
wc -l)

echo "Number of rRNA genes : $NGENES"

###############################
# rRNA featureCounts
###############################

echo
echo "Running featureCounts on rRNA genes..."

featureCounts \
-T $THREADS \
-p \
--countReadPairs \
-B \
-C \
-s $STRAND \
-t exon \
-g gene_id \
-a "$RRNA_GTF" \
-o "$OUT/rRNA_Counts.txt" \
"${BAMS[@]}"

###############################
# Calculate rRNA percentage
###############################

echo
echo "Calculating rRNA percentages..."

python3 <<EOF

import pandas as pd
import os

############################################################

gene = pd.read_csv(
    "$OUT/Cardiomyocyte_RNA_counts.txt",
    sep="\t",
    comment="#"
)

rrna = pd.read_csv(
    "$OUT/rRNA_Counts.txt",
    sep="\t",
    comment="#"
)

############################################################

# Rename columns to BAM filenames only
gene.columns = list(gene.columns[:6]) + [
    os.path.basename(c) for c in gene.columns[6:]
]

rrna.columns = list(rrna.columns[:6]) + [
    os.path.basename(c) for c in rrna.columns[6:]
]

sample_columns = gene.columns[6:]

gene_counts = gene[sample_columns].sum()

rrna_counts = rrna[sample_columns].sum()

percent = (rrna_counts / gene_counts) * 100

summary = pd.DataFrame({

    "Sample" : [x.split("/")[-1].replace(".bam","") for x in sample_columns],

    "Assigned_Gene_Counts" : gene_counts.values.astype(int),

    "Assigned_rRNA_Counts" : rrna_counts.values.astype(int),

    "rRNA_Percentage" : percent.round(3)

})

summary.to_csv(
    "$OUT/rRNA_Percentage.tsv",
    sep="\t",
    index=False
)

print(summary)

EOF

echo
echo "==========================================="
echo "Pipeline completed successfully."
date
echo "==========================================="
