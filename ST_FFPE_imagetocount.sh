#!/bin/bash

module load spaceranger/1.3.0
parentDir=/work/SCCC/s226953/spatial_experiment/05.17.2024

imageDir=$parentDir/images/
fqDir=$parentDir/ev115Lt0/


transcriptome_mm10=/work/SCCC/s226953/spatial_experiment/reference_visium_10X/refdata-gex-mm10-2020-A
transcriptome_hg38=/work/SCCC/s226953/spatial_experiment/reference_visium_10X/refdata-gex-GRCh38-2020-A
probeset_mm10=/work/SCCC/s226953/spatial_experiment/reference_visium_10X/probesets/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv
probeset_hg38=/work/SCCC/s226953/spatial_experiment/reference_visium_10X/probesets/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv


mkdir $parentDir/count
cd $parentDir/count



spaceranger count --id=ev115Lt0 \
--transcriptome=$transcriptome_hg38 \
--probe-set=$probeset_hg38 \
--fastqs=$fqDir \
--sample=ev115Lt0 \
--image=$imageDir/ev115Lto.tif \
--reorient-images \
--slide=V10S21-398 \
--area=A1 \
--localcores=20 \
--localmem=64
