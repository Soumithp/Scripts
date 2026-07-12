# Spatial transcriptomics (10x Visium)

Preprocessing scripts for 10x Visium spatial transcriptomics runs, using Space Ranger. Covers both
the FASTQ-generation step and the FFPE image-to-count step, set up for mouse (mm10) and human
(GRCh38) references with the Visium probe sets.

### Scripts

- `ST_FFPE_FASTQ.sh` — `spaceranger mkfastq` to demultiplex the BCL run into FASTQs.
- `ST_FFPE_imagetocount.sh` — `spaceranger count` on FFPE sections: takes the tissue image + FASTQs
  and produces the spot-by-gene count matrix (mm10 or GRCh38, with the matching probe set).

### Notes

- Needs the `spaceranger` module and the 10x reference transcriptomes + probe-set CSVs.
- Directory paths are from my cluster runs; update `parentDir` and the reference paths for your setup.
