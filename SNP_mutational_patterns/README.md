# Somatic/germline variants & mutational signatures

An end-to-end variant pipeline — from recalibrated BAMs through variant calling and annotation, into
either COSMIC mutational-signature fitting (somatic) or a genotype matrix (germline).
`run_snp_pipeline.sh` runs the whole thing; set `MODE=somatic` or `MODE=germline`.

### End-to-end order

1. `01_preprocess_bam.sh` — MarkDuplicates (Picard) + base quality recalibration (GATK BQSR).
2. Variant calling:
   - `02a_call_somatic_mutect2.sh` — somatic SNVs/indels, tumor vs. matched normal (Mutect2 →
     FilterMutectCalls). Somatic calls are the usual input for mutational signatures.
   - `02b_call_germline_haplotypecaller.sh` — germline joint calling (HaplotypeCaller GVCF →
     CombineGVCFs → GenotypeGVCFs). Produces the multi-sample VCF `multivcf2mat.py` expects.
3. `03_annotate_annovar.sh` — annotate with ANNOVAR (refGene, dbSNP, gnomAD, COSMIC).
4. Downstream:
   - somatic → `MutationalPatterns.R` (96-context matrix + COSMIC fit) → `normalize_mutationalpatterns.R`.
   - germline → `multivcf2mat.py` (genotype matrix); `dbSNP2mat.py` for common-variant filtering.

### Run it

```bash
# edit variables + set MODE at the top first
bash run_snp_pipeline.sh
```

### Notes

- `prepare_annovar_user.pl` and `download_prepare_annovar_user_pl.sh` are from **ANNOVAR (Kai Wang)** —
  third-party, kept with the original author credit.
- Reference genome, gnomAD, dbSNP, and COSMIC files are passed in — set them for your build (hg38, etc.).
