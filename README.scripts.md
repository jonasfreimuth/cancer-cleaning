# Information on code location and flow

## Location

* `scripts`: Main analysis runner scripts.
* `functions`: Individual functions used everywhere.
* `aux_scripts`: Additional scripts that may be necessary for getting the
  analyses ready to run.
* `notebooks`: Exploratory analysis notebooks.
* `etc`: Additional code / data that doesn't neatly fit anywhere else.

## Main workflow

1. `aux_scripts/download_dataset.sh`: Downloads & unpacks the sc-RNA-seq dataset
   into the `datasets` dir.
2. (`aux_scripts/generate_test_data.R`: If testing of analyses is intended, this
   script generates a downsampled dataset with 10% of the original cells &
   transcripts.)
3. `scripts/deconvolution_benchmark.R`: Runs the benchmark analysis of the
   current method with the specified parameters (see script for details).

## Other scripts

* `functions/benchmark_functions.R`: Main analysis functions used for carrying
  out the deconvolution and assessing the performance.
* `functions/dedupe_sigmut_mat.R`: Function copied from PiGx-SARS-CoV-2, used
  for collapsing duplicated columns in a signature matrix, mainly relevant for
  binary signature matrices.
* `functions/err_fun_common.R`: Code for performing tasks for all the error
  functions used in the project.
* `functions/norm_functions.R`: Functions for normalizing the count data used.
* `functions/rmse.R`: Function for computing the root mean squared error of
  two vectors.
* `functions/util_functions.R`: Various utility functions used by the analysis &
  auxiliary scripts.

## Unused scripts

None at the moment. Generally, unused scripts are removed and specific analysis
artefacts should always be associated with the code that produced them in some
way.
