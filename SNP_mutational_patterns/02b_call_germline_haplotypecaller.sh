#!/bin/bash
# Germline variant calling with GATK HaplotypeCaller (GVCF) + joint genotyping.
# Produces the multi-sample VCF that multivcf2mat.py expects (GenotypeGVCFs output).
module load gatk/4.2.6.1

ref=/path/to/genome.fa
bamDir=/path/to/bam_recal
outDir=/path/to/germline_vcf
mkdir -p "$outDir"

# per-sample GVCFs
for bam in "$bamDir"/*.recal.bam; do
    sample=$(basename "$bam" .recal.bam)
    echo "HaplotypeCaller: $sample"
    gatk HaplotypeCaller \
        -R "$ref" -I "$bam" \
        -ERC GVCF \
        -O "$outDir/${sample}.g.vcf.gz"
done

# combine + joint genotype across the cohort
gvcfArgs=()
for g in "$outDir"/*.g.vcf.gz; do gvcfArgs+=(-V "$g"); done

gatk CombineGVCFs -R "$ref" "${gvcfArgs[@]}" -O "$outDir/cohort.g.vcf.gz"
gatk GenotypeGVCFs -R "$ref" -V "$outDir/cohort.g.vcf.gz" -O "$outDir/cohort.joint.vcf.gz"
echo "joint-genotyped VCF -> $outDir/cohort.joint.vcf.gz"
