#!/bin/bash
# BAM preprocessing before variant calling: mark duplicates + base quality recalibration.
module load picard/2.26
module load gatk/4.2.6.1
module load samtools/1.15

ref=/path/to/genome.fa            # indexed (.fai) with a .dict alongside
dbsnp=/path/to/dbsnp.vcf.gz       # known sites for BQSR
bamDir=/path/to/bam
outDir=/path/to/bam_recal
mkdir -p "$outDir"

for bam in "$bamDir"/*.bam; do
    sample=$(basename "$bam" .bam)
    echo "preprocessing $sample"

    picard MarkDuplicates \
        I="$bam" \
        O="$outDir/${sample}.dedup.bam" \
        M="$outDir/${sample}.dup_metrics.txt"
    samtools index "$outDir/${sample}.dedup.bam"

    gatk BaseRecalibrator \
        -I "$outDir/${sample}.dedup.bam" \
        -R "$ref" --known-sites "$dbsnp" \
        -O "$outDir/${sample}.recal.table"

    gatk ApplyBQSR \
        -I "$outDir/${sample}.dedup.bam" \
        -R "$ref" --bqsr-recal-file "$outDir/${sample}.recal.table" \
        -O "$outDir/${sample}.recal.bam"
done
echo "recalibrated BAMs -> $outDir"
