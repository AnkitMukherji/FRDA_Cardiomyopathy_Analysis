#!/bin/bash

# 'adapterremoval', 'afterqc', 'anglerfish', 'ataqv'
# 'bakta', 'bamdst', 'bamtools', 'bases2fastq', 'bbduk', 'bbmap', 'bcftools', 'bcl2fastq', 'bclconvert', 'biobambam2', 'biobloomtools', 'biscuit', 'bismark', 'bowtie1', 'bowtie2', 'busco', 'bustools'
# 'ccs', 'cellranger', 'cellranger_arc', 'cells2stats', 'checkatlas', 'checkm', 'checkm2', 'checkqc', 'clipandmerge', 'clusterflow', 'conpair', 'custom_content', 'cutadapt'
# 'damageprofiler', 'deacon', 'dedup', 'deeptools', 'diamond', 'disambiguate', 'dragen', 'dragen_fastqc'
# 'eigenstratdatabasetools'
# 'fastp', 'fastq_screen', 'fastqc', 'fastqe', 'featurecounts', 'fgbio', 'filtlong', 'flash', 'flexbar', 'freyja'
# 'ganon', 'gatk', 'gffcompare', 'glimpse', 'goleft_indexcov', 'gopeaks', 'gtdbtk'
# 'haplocheck', 'happy', 'hicexplorer', 'hicpro', 'hicstuff', 'hicup', 'hifi_trimmer', 'hifiasm', 'hisat2', 'homer', 'hops', 'hostile', 'htseq', 'humid'
# 'interop', 'isoseq', 'ivar'
# 'jcvi', 'jellyfish'
# 'kaiju', 'kallisto', 'kat', 'kraken'
# 'leehom', 'librarian', 'lima', 'longranger'
# 'macs2', 'malt', 'mapdamage', 'megahit', 'metaphlan', 'methurator', 'methylqa', 'mgikit', 'minionqc', 'mirtop', 'mirtrace', 'mosaicatcher', 'mosdepth', 'motus', 'mtnucratio', 'multivcfanalyzer'
# 'nanoq', 'nanostat', 'nextclade', 'ngsbits', 'ngsderive', 'nonpareil'
# 'odgi', 'optitype'
# 'pairtools', 'pangolin', 'pbmarkdup', 'peddy', 'percolator', 'phantompeakqualtools', 'picard', 'porechop', 'preseq', 'prinseqplusplus', 'prokka', 'purple', 'pychopper', 'pycoqc'
# 'qc3C', 'qorts', 'qualimap', 'quast'
# 'ribotish', 'ribowaltz', 'riker', 'rna_seqc', 'rockhopper', 'rsem', 'rseqc'
# 'salmon', 'sambamba', 'samblaster', 'samtools', 'sargasso', 'seqera_cli', 'seqfu', 'seqkit', 'sequali', 'seqwho', 'seqyclean', 'sexdeterrmine', 'sickle', 'sincei', 'skewer', 'slamdunk', 'snippy', 'snpeff', 'snpsplit', 'somalier', 'sompy', 'sortmerna', 'sourmash', 'spaceranger', 'stacks', 'star', 'supernova', 'sylphtax'
# 'telseq', 'theta2', 'tophat', 'trim_galore', 'trimmomatic', 'truvari'
# 'umicollapse', 'umitools'
# 'varscan2', 'vcftools', 'vep', 'verifybamid', 'vg'
# 'whatshap'
# 'xengsort', 'xenium', 'xenome'

OUT="/mnt/faruq2/lab_users/ankit/rna_seq_asangla/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_MultiQC"

mkdir -p "$OUT"

multiqc \
    --module dragen \
    --module dragen_fastqc \
    --module fastqc \
    --module fastp \
    --module rseqc \
    --module featurecounts \
    --module qualimap \
    -o "$OUT" \
    --filename RNA_MultiQC_Report_Cardiomyocyte_Day52 \
    --fn_as_s_name \
    -f \
    /mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Fastq_Reverse \
    /mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Fastq_Reverse_Trimmed \
    /mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Dragen_Mapping_Fastq_Reverse_Trimmed \
    /mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Strandedness \
    /mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Feature_Counts \
    /mnt/faruq2/lab_users/ankit/rna_seq_cardio/rna_seq_asangla_cardio_12_06_26_lane_001_day_52/RNA_Qualimap_Results
