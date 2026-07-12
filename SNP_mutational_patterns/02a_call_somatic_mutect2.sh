#!/bin/bash
# Somatic variant calling with GATK Mutect2 (tumor vs matched normal).
# Somatic SNVs/indels are the usual input for tumor mutational-signature analysis.
module load gatk/4.2.6.1

ref=/path/to/genome.fa
gnomad=/path/to/af-only-gnomad.vcf.gz     # germline resource
pairsFile=/path/to/tumor_normal_pairs.tsv # columns: sample  tumor.bam  normal.bam  normalSampleName
outDir=/path/to/somatic_vcf
mkdir -p "$outDir"

# one tumor/normal pair per line
while IFS=$'\t' read -r sample tumorBam normalBam normalName; do
    [ -z "$sample" ] && continue
    echo "calling somatic variants: $sample"
    gatk Mutect2 \
        -R "$ref" \
        -I "$tumorBam" -I "$normalBam" -normal "$normalName" \
        --germline-resource "$gnomad" \
        -O "$outDir/${sample}.unfiltered.vcf.gz"

    gatk FilterMutectCalls \
        -R "$ref" \
        -V "$outDir/${sample}.unfiltered.vcf.gz" \
        -O "$outDir/${sample}.filtered.vcf.gz"
done < "$pairsFile"
echo "somatic VCFs -> $outDir"
