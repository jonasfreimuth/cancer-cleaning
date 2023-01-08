# TODO Overhaul the output path handling / creating
## ----setup--------------------------------------------------------------------
library("magrittr")
library("here")

here::i_am("scripts/runner.R")

source(here("functions/util_functions.R"))

source(here("R6/ScRnaExperiment.R"))
source(here("R6/ReferenceParams.R"))
source(here("R6/Reference.R"))
source(here("R6/PseudobulkParams.R"))
source(here("R6/Pseudobulk.R"))
source(here("R6/DeconvParams.R"))
source(here("R6/Deconvolution.R"))
source(here("R6/DeconvSumParams.R"))
source(here("R6/DeconvolutionSummary.R"))

## ----parameters --------------------------------------------------------------
# TODO Clean up args after normalization fixing.
cmd_args <- commandArgs(trailingOnly = TRUE)

# names must match script params
arg_names <- c(
  # From where ground truth data should be taken. The root of that directory
  # needs to contain the files
  # * count_matrix_sparse.mtx
  # * count_matrix_genes.tsv
  # * count_matrix_barcodes.tsv
  # * "metadata.csv"
  "data_path",

  # Possible types:
  # * binary: Produces a binary sigmat based on the threshold
  # * deseq2: Selects the top n differentially expressed genes across the whole
  #   count_mat, as determined by DESeq2. n is based on the threshold.
  "sigmat_type",

  # Arguments for the sequence exploring the threshold parameter for the
  # sigmat type. Everything except the power should be a fraction.
  # Power gives the power law by which elements of the sequence are spaced.
  "thresh_start",
  "thresh_stop",
  "thresh_step",
  "thresh_power",

  # Number of pseudobulks
  "n_pseudobulk",

  # Fraction of count_mat cells sampled per pseudobulk
  "pseudobulk_cell_frac",

  # Possible types:
  # * norm: Plain normalization, each per cell transcript count divided by total
  #   per transcript counts.
  # * lognorm: Same as norm, but the base 10 log is taken for each individual
  #   count.
  # * (quantile: Quantile normalization of transcript counts across cells.
  #    Will probably lead to weird results.)
  "normalization_type",

  # Whether normalization should be just once applied to the count mat, or to
  # pseudobulk and signature matrix independently. The independent option is
  # recommended, this option will probably be removed soon.
  "normalize_independently",

  # The pseudobulk deconvolution method to be used. See deconvR for details.
  "deconv_method",

  # Seed used for random sampling of pseudobulks.
  "seed"
)

default_args <- c(
  data_path = "datasets/Wu_etal_downsampled_test/",
  thresh_start = 0.5,
  thresh_stop = 0.5,
  thresh_step = 0.1,
  thresh_power = 1,
  sigmat_type = "deseq2",
  n_pseudobulk = "10",
  pseudobulk_cell_frac = "0.2",
  normalization_type = "lognorm",
  normalize_independently = "TRUE",

  # TODO Investigate why QP fails.
  deconv_method = "nnls",
  seed = "123"
)

# Priority:
#   1. commandArgs (If provided correctly, i.e. the right length)
#   2. script_args in the workspace (So that manual specification of non-default
#      args can persist)
#   3. default_args (As fallback.)
if (length(cmd_args) == length(arg_names)) {
  cat("\nUsing provided args...\n")
  script_args <- cmd_args
  names(script_args) <- arg_names
} else if (exists("script_args")) {
  cat("\nUsing existing script_args variable...\n")
} else {
  cat("\nUsing default args...\n")
  script_args <- default_args
}

# Change to list for nicer access to params. params will also include some
# parameters that are to be tweaked just from within the script.
params <- as.list(script_args)

# Ensure correct types
# TODO Make this vectorized somehow, maybe give the arg_names vector the name
# of the coercion function to be applied or something.
# TODO Deal with params that might evaluate to something else than an atomic
# vector.
params$normalize_independently %<>%
  as.logical()
params$n_pseudobulk %<>%
  as.numeric()
params$pseudobulk_cell_frac %<>%
  as.numeric()
params$seed %<>%
  as.numeric()
params$thresh_start %<>%
  as.numeric()
params$thresh_stop %<>%
  as.numeric()
params$thresh_step %<>%
  as.numeric()
params$thresh_power %<>%
  as.numeric()
params$sigmat_type %<>%
  tolower()

params <- lapply(
  params,
  function(param) {
    if (length(param) > 0) {
      if (!is.null(param) & !is.na(param)) {
        return(param)
      }
    }
  }
)

params$norm_scale <- 1

params$base_width <- 3
params$base_height <- 2

params$facet_height <- 2
params$facet_width <- 5.5

param_sep <- "_"
pair_sep <- "-"

parameter_string <- paste(
  paste0("data_path", pair_sep, basename(params$data_path)),
  paste0("seed", pair_sep, params$seed),
  paste0("method", pair_sep, params$deconv_method),
  paste0("indepnorm", pair_sep, params$normalize_independently),
  paste0("normtype", pair_sep, params$normalization_type),
  paste0("npseudobulk", pair_sep, params$n_pseudobulk),
  paste0("sigmat_type", pair_sep, params$sigmat_type),
  paste0("sizefrac", pair_sep, params$pseudobulk_cell_frac),
  sep = param_sep
)

run_path <- here(
  "output",
  paste(
    format(Sys.time(), "%Y%m%d-%H%M%S"),
    parameter_string,
    sep = param_sep
  )
)

dir.create(run_path, recursive = TRUE, showWarnings = FALSE)

# Create executable record of the args used in this script
sink(file = here(run_path, "args.R"))
cat("script_args <- c(\n\t")
script_args %>%
  {
    paste0(names(.), " = \"", ., "\"")
  } %>%
  cat(sep = ",\n\t")
cat(")")
sink()

param_txt <- params %>%
  unlist() %>%
  {
    paste(names(.), ., sep = ": ")
  } %>%
  paste(collapse = "\n") %>%
  paste0("\n")

# Display params
cat(param_txt)

# Save final params
cat(param_txt, file = here(run_path, "params.txt"))

set.seed(params$seed)

if (params$normalize_independently) {
  super_norm_type <- NULL
  sub_norm_type <- params$normalization_type
} else {
  super_norm_type <- params$normalization_type
  sub_norm_type <- NULL
}


# Bundle the code used for the analysis with the results. Seems alright as long
# as it is not too much code.
file.copy(
  here("manifest.scm"),
  here(run_path, "manifest.scm")
) %>%
  invisible()

file.copy(
  here("scripts/deconvolution_benchmark.R"),
  here(run_path, "deconvolution_benchmark.R")
) %>%
  invisible()

fun_target_path <- here(run_path, "functions")
dir.create(fun_target_path, recursive = TRUE, showWarnings = FALSE)

lapply(
  dir("functions", full.names = TRUE),
  function(file_path) {
    file.copy(
      file_path,
      here(fun_target_path, basename(file_path))
    )
  }
) %>%
  invisible()

fun_target_path <- here(run_path, "R6")
dir.create(fun_target_path, recursive = TRUE, showWarnings = FALSE)

lapply(
  dir("R6", full.names = TRUE),
  function(file_path) {
    file.copy(
      file_path,
      here(fun_target_path, basename(file_path))
    )
  }
) %>%
  invisible()


## ----params_to_explore -------------------------------------------------------
thresh_list <- seq_power(
  start = params$thresh_start,
  stop = params$thresh_stop,
  step = params$thresh_step,
  power = params$thresh_power
) %>%
  set_names(.) %>%
  as.list()

cancer_cols_list <- list(
  "Cancer split" = c("Cancer Epithelial"),
  "Cancer not split" = NULL
)


## ----data_loading-------------------------------------------------------------
experiment <- ScRnaExperiment$new(
  count_mat_file = here(
    params$data_path, "count_matrix_sparse.mtx"
  ),
  rowname_file = here(
    params$data_path, "count_matrix_genes.tsv"
  ),
  colname_file = here(
    params$data_path, "count_matrix_barcodes.tsv"
  ),
  meta_file = here(
    params$data_path, "metadata.csv"
  )
)


## ----signature_matrix_generation----------------------------------------------
deconv_ref_list <- lapply(
  thresh_list,
  function(thresh) {
    lapply(
      cancer_cols_list,
      function(cancer_cols) {
        experiment$create_reference(
          params = ReferenceParams$new(
            threshold = thresh,
            metric = params$sigmat_type,
            cancer_celltypes = cancer_cols,
            norm_type = params$normalization_type,
            norm_scale_factor = params$norm_scale
          )
        )
      }
    )
  }
) %>%
  unlist() %>%
  as.list()

is_null_reference <- lapply(deconv_ref_list, is.null) %>%
  unlist()

is_empty_reference <- lapply(
  deconv_ref_list,
  function(ref) {
    ref$n_transcripts <= 0
  }
) %>%
  unlist()

deconv_ref_list <- deconv_ref_list %>%
  magrittr::extract(!is_null_reference & !is_empty_reference)

# Generate heatmaps
heatmap_path <- here(run_path, "plots/heatmaps")
dir.create(heatmap_path, recursive = TRUE, showWarnings = FALSE)

lapply(
  deconv_ref_list,
  function(reference) {
    reference$print_heatmap(heatmap_path)
  }
)

## ----pseudobulk_generation----------------------------------------------------
n_bulk_cells <- params$pseudobulk_cell_frac * experiment$n_cells

i <- 0

pseudobulk_list <-
  # predraw which cells are used in each pseudobulk
  lapply(
    rep(list(seq_len(experiment$n_cells)), params$n_pseudobulk),
    sample,
    n_bulk_cells
  ) %>%
  # actually draw pseudobulks
  lapply(
    function(idx_vec) {
      i <<- i + 1
      lapply(
        cancer_cols_list,
        function(cancer_cols) {
          experiment$create_pseudobulk(
            PseudobulkParams$new(
              id = i,
              cell_indices = idx_vec,
              cancer_celltypes = cancer_cols,
              norm_type = params$normalization_type,
              norm_scale_factor = params$norm_scale
            )
          )
        }
      )
    }
  ) %>%
  unlist() %>%
  as.list()

rm(i)

## ----deconvolution------------------------------------------------------------
deconv_summary <- expand_grid(
  reference = deconv_ref_list,
  pseudobulk = pseudobulk_list
) %>%
  apply(
    1,
    function(deconv_pair) {
      ref <- deconv_pair$reference
      pbulk <- deconv_pair$pseudobulk
      # FIXME Make pseudobulk and reference matching more elegant.
      if (isTRUE(all.equal(
        sort(ref$params$cancer_celltypes),
        sort(pbulk$params$cancer_celltypes)
      ))) {
        return(
          Deconvolution$new(
            reference = deconv_pair$reference,
            pseudobulk = deconv_pair$pseudobulk,
            params = DeconvParams$new(
              deconvolution_method = params$deconv_method
            )
          )
        )
      }
    }
  ) %>%
  magrittr::extract(!unlist(lapply(., is.null))) %>%
  {
    DeconvolutionSummary$new(
      deconvolution_list = .,
      params = DeconvSumParams$new()
    )
  }

plot_path <- here(run_path, "plots")
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

deconv_summary$print_rmse_plot(dir = plot_path)