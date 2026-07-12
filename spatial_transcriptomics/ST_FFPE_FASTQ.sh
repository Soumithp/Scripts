#!/bin/bash

##################################################################################################
module load spaceranger/1.3.0
parentDir=/project/SCCC/Hoshida_lab/soumith/ST/
cd $parentDir

spaceranger mkfastq --id=fastq \
--run=$parentDir/L002_RGL25a_bcl \
--csv=$parentDir/spaceranger-simple.csv
