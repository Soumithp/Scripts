#!/bin/bash
# STAR alignment: FASTQ -> sorted BAM (primary aligner for the RNA-seq pipeline)
# Builds the genome index once, then aligns each paired-end sample.
module load star/2.7.9a
module load samtools/1.15

threads=8
genomeDir=/path/to/star_index
genomeFasta=/path/to/genome.fa          # e.g. GRCh38.primary_assembly.genome.fa
gtf=/path/to/annotation.gtf             # e.g. gencode.v44.annotation.gtf
fastqDir=/path/to/fastq
outDir=/path/to/aligned
mkdir -p "$outDir"

# --- build index once (skip if $genomeDir already built) ---
if [ ! -f "$genomeDir/SAindex" ]; then
    mkdir -p "$genomeDir"
    STAR --runMode genomeGenerate \
         --genomeDir "$genomeDir" \
         --genomeFastaFiles "$genomeFasta" \
         --sjdbGTFfile "$gtf" \
         --sjdbOverhang 100 \
         --runThreadN "$threads"
fi

# --- align each sample (expects ${sample}_R1.fastq.gz / _R2.fastq.gz) ---
for r1 in "$fastqDir"/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="$fastqDir/${sample}_R2.fastq.gz"
    echo "aligning $sample"
    STAR --runMode alignReads \
         --genomeDir "$genomeDir" \
         --readFilesIn "$r1" "$r2" \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         --runThreadN "$threads" \
         --outFileNamePrefix "$outDir/${sample}_"
    samtools index "$outDir/${sample}_Aligned.sortedByCoord.out.bam"
done
echo "STAR alignment done -> $outDir"
