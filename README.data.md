# Information on data used in this project

## Raw data

Currently, this project exclusively uses sc-RNA-seq data from Wu, S.Z.,
Al-Eryani, G., Roden, D.L. et al. A single-cell and spatially resolved atlas of
human breast cancers. Nat Genet 53, 1334–1347 (2021).
<https://doi.org/10.1038/s41588-021-00911-1>, downloaded from
<ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz>
via `aux_scripts/download_dataset.sh`. As this data is freely available and only
2.3 GB large unpacked, the project is set up for the data to be located in
the `datasets` dir inside.

## Processed data

No intermediate files are currently produced.

## Tertiary data

All analysis output is saved to the `output` dir, with each run having its own
timestamped dir. The latter generally contains:

* `functions`: A dir with a copy of all the function files used in generating
  the results.
* `plots`: A dir containing all the output plots, currently those are
  * `pred_prop_expr_all_plot.png`: Plot of transcript proportions predicted by
    deconvolution versus true cancer expression, for all transcripts.
  * `pred_prop_expr_marker_plot.png`: Plot of transcript proportions predicted
    by deconvolution versus true cancer expression, for marker transcripts only.
  * `resid_expr_all_plot.png`: Plot of deconvolution transcript residuals versus
    true cancer expression, for all transcripts.
  * `resid_expr_marker_plot.png`: Plot of deconvolution transcript residuals
    versus true cancer expression, for marker transcripts only.
  * `rmse_plot.png`: Plot of per celltype root mean squared errors.
* `args.R`: A sourceable R file containing the arguments to the main script used
  to run the anaylsis.
* `cancer_comparison_summary.csv`: A summary table of used parameters, and
  performance metrics for the analysis.
* `deconvolution_benchmark.R`: A copy of the main analysis script.
* `params.txt`: A plain text file giving all of the parameters used in the
  script (i.e. args & some parameters defined in the script itself).
