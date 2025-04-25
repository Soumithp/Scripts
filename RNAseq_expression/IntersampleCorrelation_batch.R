setwd("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/03.15.2024_RGL")
source("IntersampleCorrelation.R")

IntersampleCorrelation(
    "RGL_2_RLE.gct",
    reference.method="geo.mean", # or "mean", "geo.mean"
    add.one="T",
    output.name="IntCor_RGL_2"
    )


