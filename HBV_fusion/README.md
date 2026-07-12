# HBV fusion / viral integration detection

Scripts I used to find hepatitis B virus (HBV) integration and viral–host fusion events in
RNA-seq data, mainly in the context of hepatocellular carcinoma. The detection step is built on
Kraken (for viral read filtering) and ViFi / FastViFi (for integration calling), wrapped so it can
be run per sample and then summarized across a cohort.

### Rough order

1. `build_ref_idx.sh` — set up the ViFi/AmpliconArchitect data repo and reference indices.
2. `build_kraken_hbv.sh` — build the Kraken database used to pull out viral reads.
3. `run_kraken_vifi_pipeline.py` — the main per-sample driver. Kraken-filters reads, then runs ViFi
   to call integration events. Supports several viruses (hbv, hcv, hpv, ebv) and a
   sample-level vs. sensitive mode.
4. `hbvfusion_sum.py` — turns the raw integration calls into tidy tables: human gene/exon on one
   side, viral gene on the other, plus a merged junction matrix.
5. `cluster_trans_new.py` — collapses/clusters fusion transcripts that map to the same event.
6. `norm_fusion_readc.R` — normalizes fusion read counts across samples.

### Notes

- `download_prepare_annovar_user_pl.sh` fetches a helper from **ANNOVAR** (Kai Wang) — third-party,
  used for annotation. Kept as-is with its original attribution.
- Paths in the scripts point to where I ran these on the cluster; change them to your own before running.
