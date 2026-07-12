#!/bin/bash
# Annotate a VCF with ANNOVAR (gene, region, and filter-based databases).
# ANNOVAR is third-party (Kai Wang); see prepare_annovar_user.pl in this folder.
module load annovar

humandb=/path/to/annovar/humandb        # downloaded ANNOVAR databases
buildver=hg38
inVcf=${1:?usage: 03_annotate_annovar.sh <input.vcf.gz> <out_prefix>}
outPrefix=${2:?usage: 03_annotate_annovar.sh <input.vcf.gz> <out_prefix>}

table_annovar.pl "$inVcf" "$humandb" \
    -buildver "$buildver" \
    -out "$outPrefix" \
    -remove \
    -protocol refGene,cytoBand,avsnp150,gnomad211_exome,cosmic70 \
    -operation g,r,f,f,f \
    -nastring . \
    -vcfinput
echo "annotated -> ${outPrefix}.${buildver}_multianno.txt"
