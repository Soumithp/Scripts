#!/bin/bash
# HISAT2 alignment: FASTQ -> sorted BAM (lighter alternative to STAR)
module load hisat2/2.2.1
module load samtools/1.15

threads=8
hisat2Index=/path/to/hisat2_index/genome   # prefix built with hisat2-build
fastqDir=/path/to/fastq
outDir=/path/to/aligned_hisat2
mkdir -p "$outDir"

for r1 in "$fastqDir"/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="$fastqDir/${sample}_R2.fastq.gz"
    echo "aligning $sample"
    hisat2 -p "$threads" -x "$hisat2Index" -1 "$r1" -2 "$r2" \
        | samtools sort -@ "$threads" -o "$outDir/${sample}.sorted.bam"
    samtools index "$outDir/${sample}.sorted.bam"
done
echo "HISAT2 alignment done -> $outDir"
