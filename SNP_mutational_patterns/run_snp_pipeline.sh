#!/bin/bash
# End-to-end SNP / mutational-signature pipeline.
# preprocess BAM -> variant calling -> annotation -> (somatic) mutational signatures
#                                                  \-> (germline) genotype matrix
# Set MODE to "somatic" or "germline", fill the variables, then: bash run_snp_pipeline.sh
set -euo pipefail

MODE="somatic"                     # "somatic" (Mutect2) or "germline" (HaplotypeCaller)
ref=/path/to/genome.fa
cosmicSig=/path/to/COSMIC_signatures.txt
refCode="hg38"
workDir=/path/to/work
scriptDir=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$workDir"

# ---- 1. BAM preprocessing (dedup + BQSR) ----
bash "$scriptDir/01_preprocess_bam.sh"

if [ "$MODE" = "somatic" ]; then
    # ---- 2. somatic calling (Mutect2) ----
    bash "$scriptDir/02a_call_somatic_mutect2.sh"

    # ---- 3. annotate + 4. per-sample mutational signatures ----
    for vcf in "$workDir"/somatic_vcf/*.filtered.vcf.gz; do
        sample=$(basename "$vcf" .filtered.vcf.gz)
        bash "$scriptDir/03_annotate_annovar.sh" "$vcf" "$workDir/annot/${sample}"
        Rscript --vanilla "$scriptDir/MutationalPatterns.R" \
            "$vcf" "$sample" "$refCode" "$cosmicSig"
    done
    Rscript --vanilla "$scriptDir/normalize_mutationalpatterns.R"

else
    # ---- 2. germline joint calling (HaplotypeCaller -> GenotypeGVCFs) ----
    bash "$scriptDir/02b_call_germline_haplotypecaller.sh"
    jointVcf="$workDir/germline_vcf/cohort.joint.vcf.gz"

    # ---- 3. annotate + 4. genotype matrix ----
    bash "$scriptDir/03_annotate_annovar.sh" "$jointVcf" "$workDir/annot/cohort"
    python "$scriptDir/multivcf2mat.py" --vcf "$jointVcf" --out "$workDir/genotype_matrix.txt" --DP 10
fi

echo "SNP pipeline ($MODE) finished"
