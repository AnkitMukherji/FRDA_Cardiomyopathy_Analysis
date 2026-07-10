#!/bin/bash

INPUT="/mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Qualimap_Results"
OUTPUT="/mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Qualimap_Results"
TSV_OUT="$OUTPUT/RNA_Asangla_Day52_Qualimap_Summary.tsv"

printf "Sample\tExonic_Millions\tExonic_%%\tIntronic_Millions\tIntronic_%%\tIntergenic_Millions\tIntergenic_%%\n" > "$TSV_OUT"

for dir in "$INPUT"/*/; 
do
    sample="${dir%/}"
    sample="${sample##*/}"
    file=$(find "$dir" -maxdepth 1 -name "*.txt" | head -n 1)
    
    echo "Found Sample Name : $sample"
    if [ -n "$file" ]; then
        echo "Found File Path   : $file"
    fi
    
    if [ -f "$file" ]; then
        metrics=$(awk '
            /exonic|intronic|intergenic/ {
                gsub(/,/, "", $3) # Removes , to convert into numeric
                gsub(/[()%]/, "", $4) # Extracts the number removing percentage
                val[$1] = sprintf("%.3f\t%s", $3/1000000, $4) # %s is for placeholder text
            }
            END { print val["exonic"] "\t" val["intronic"] "\t" val["intergenic"] }
        ' "$file")
        
        printf "%s\t%s\n" "$sample" "$metrics" >> "$TSV_OUT"
    fi
done

