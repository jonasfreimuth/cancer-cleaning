# TODO Overhaul the output path handling / creating
## ----setup--------------------------------------------------------------------
library("magrittr")
library("here")

here::i_am("scripts/deconvolution_benchmark.R")

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
  "Cancer split" = c("Cancer Celltypes"),
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
      i <- i + 1
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

# ## ----plot_deconv_err----------------------------------------------------------
# sample_sum_df <- split_res_list %>%
#   dfextract2("deconv_sum", "sigmat_thresh", "split")
#
# rename_vec <- sample_sum_df %>%
#   head(n = 1) %>%
#   select(matches(c("rmse", "cancer_expr_v_*"))) %>%
#   names()
#
# rename_targets <- rename_vec %>%
#   str_replace_all(pattern = "cancer_expr", replacement = "cexpr") %>%
#   str_replace_all(pattern = "rmse", replacement = "bulk_rmse")
#
# rename_vec <- rename_vec %>%
#   set_names(rename_targets)
#
# parameter_sum_df <- sample_sum_df %>%
#   group_by(split, sigmat_thresh) %>%
#   rename(!!!rename_vec) %>%
#   summarize(
#     across(
#       matches(c("bulk_rmse", "cexpr_v_*")),
#       list("mean" = mean),
#       .names = "{.fn}_{.col}"
#     ),
#     .groups = "drop_last"
#   )
#
#
# celltype_rmse_df <- split_res_list %>%
#   dfextract2("deconv_res", "sigmat_thresh", "split")
#
# celltype_df_sum <- celltype_rmse_df %>%
#   group_by(split, sigmat_thresh, celltype) %>%
#   summarise(
#     # per celltype rmse, across all samples
#     celltype_rmse = rmse(prop_true, prop_deconv),
#     .groups = "drop_last"
#   ) %>%
#   left_join(parameter_sum_df, by = c("split", "sigmat_thresh")) %>%
#   mutate(sigmat_thresh = as.numeric(sigmat_thresh))
#
#
# plot_err <- ggplot(
#   celltype_df_sum,
#   aes(celltype, celltype_rmse)
# ) +
#   geom_col(alpha = 0.5, position = position_identity()) +
#   geom_text(
#     aes(
#       x = length(unique(celltype)) / 2,
#       y = max(celltype_rmse) * 1.1,
#       label = paste0(
#         "Mean bulk RMSE: ", round(mean_bulk_rmse, 3)
#       )
#     ),
#     vjust = 1
#   ) +
#   labs(
#     x = "Cell types",
#     y = "Per celltype RMSE"
#   ) +
#   facet_grid(
#     cols = vars(split),
#     rows = vars(sigmat_thresh)
#   ) +
#   theme_benchmark() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
# dir.create(here(run_path, "plots"), recursive = TRUE, showWarnings = FALSE)
#
# ggsave(here(run_path, "plots", "rmse_plot.png"), plot_err,
#   height = params$base_height + (
#     length(thresh_seq) * params$facet_height
#   ),
#   width = params$base_width + (
#     length(split_vals) * params$facet_width
#   )
# )
#
# # Resid vs Expr Plot ------------------------------------------------------
# resid_expr_marker_df <- split_res_list %>%
#   dfextract2("resid_expr_marker_df", "sigmat_thresh", "split") %>%
#   left_join(parameter_sum_df, by = c("split", "sigmat_thresh"))
#
# plot_resid_expr_marker <- ggplot(
#   resid_expr_marker_df,
#   aes(cancer_expr, resid, col = sample)
# ) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_text(
#     aes(
#       x = mean(c(min(cancer_expr), max(cancer_expr))),
#       y = max(resid) * 1.1,
#       label = paste0(
#         "Mean r: ", round(mean_cexpr_v_resid_marker, 3)
#       )
#     ),
#     vjust = 1,
#     col = "gray33"
#   ) +
#   labs(
#     x = "Cancer expression",
#     y = "Transcript residual"
#   ) +
#   facet_wrap(~sample) +
#   facet_grid(
#     cols = vars(split),
#     rows = vars(sigmat_thresh),
#     scales = "free"
#   ) +
#   theme_benchmark() +
#   theme(legend.position = "none")
#
# dir.create(here(run_path, "plots"), recursive = TRUE, showWarnings = FALSE)
#
# ggsave(here(run_path, "plots", "resid_expr_marker_plot.png"),
#   plot_resid_expr_marker,
#   height = params$base_height + (
#     length(thresh_seq) * params$facet_height
#   ),
#   width = params$base_width + (
#     length(split_vals) * params$facet_width
#   )
# )
#
# resid_expr_all_df <- split_res_list %>%
#   dfextract2("resid_expr_all_df", "sigmat_thresh", "split") %>%
#   left_join(parameter_sum_df, by = c("split", "sigmat_thresh"))
#
# plot_resid_expr_all <- ggplot(
#   resid_expr_all_df,
#   aes(cancer_expr, resid, col = sample)
# ) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_text(
#     aes(
#       x = mean(c(min(cancer_expr), max(cancer_expr))),
#       y = max(resid) * 1.1,
#       label = paste0(
#         "Mean r: ", round(mean_cexpr_v_resid_all, 3)
#       )
#     ),
#     vjust = 1,
#     col = "gray33"
#   ) +
#   labs(
#     x = "Cancer expression",
#     y = "Transcript residual"
#   ) +
#   facet_wrap(~sample) +
#   facet_grid(
#     cols = vars(split),
#     rows = vars(sigmat_thresh),
#     scales = "free"
#   ) +
#   theme_benchmark() +
#   theme(legend.position = "none")
#
# dir.create(here(run_path, "plots"), recursive = TRUE, showWarnings = FALSE)
#
# ggsave(here(run_path, "plots", "resid_expr_all_plot.png"), plot_resid_expr_all,
#   height = params$base_height + (
#     length(thresh_seq) * params$facet_height
#   ),
#   width = params$base_width + (
#     length(split_vals) * params$facet_width
#   )
# )
#
#
# # Pred Trans Abund vs Expr Plot -------------------------------------------
# plot_pred_prop_expr_marker <- ggplot(
#   resid_expr_marker_df,
#   aes(cancer_expr, deconv_pred, col = sample)
# ) +
#   geom_point() +
#   geom_text(
#     aes(
#       x = mean(c(min(cancer_expr), max(cancer_expr))),
#       y = max(deconv_pred) * 1.2,
#       label = paste0(
#         "Mean r: ", round(mean_cexpr_v_deconv_pred_marker, 3)
#       )
#     ),
#     vjust = 1,
#     col = "gray33"
#   ) +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(
#     x = "Cancer expression",
#     y = "Deconvolution model prediction"
#   ) +
#   facet_wrap(~sample) +
#   facet_grid(
#     cols = vars(split),
#     rows = vars(sigmat_thresh),
#     scales = "free"
#   ) +
#   theme_benchmark() +
#   theme(legend.position = "none")
#
# ggsave(here(run_path, "plots", "pred_prop_expr_marker_plot.png"),
#   plot_pred_prop_expr_marker,
#   height = params$base_height + (
#     length(thresh_seq) * params$facet_height
#   ),
#   width = params$base_width + (
#     length(split_vals) * params$facet_width
#   )
# )
#
#
# plot_pred_prop_expr_all <- ggplot(
#   resid_expr_all_df,
#   aes(cancer_expr, deconv_pred, col = sample)
# ) +
#   geom_point() +
#   geom_text(
#     aes(
#       x = mean(c(min(cancer_expr), max(cancer_expr))),
#       y = max(deconv_pred) * 1.2,
#       label = paste0(
#         "Mean r: ", round(mean_cexpr_v_deconv_pred_all, 3)
#       )
#     ),
#     vjust = 1,
#     col = "gray33"
#   ) +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(
#     x = "Cancer expression",
#     y = "Deconvolution model prediction"
#   ) +
#   facet_wrap(~sample) +
#   facet_grid(
#     cols = vars(split),
#     rows = vars(sigmat_thresh),
#     scales = "free"
#   ) +
#   theme_benchmark() +
#   theme(legend.position = "none")
#
# ggsave(here(run_path, "plots", "pred_prop_expr_all_plot.png"),
#   plot_pred_prop_expr_all,
#   height = params$base_height + (
#     length(thresh_seq) * params$facet_height
#   ),
#   width = params$base_width + (
#     length(split_vals) * params$facet_width
#   )
# )
#
# # ----Data saving---------------------------------------------------------------
# print(parameter_sum_df)
#
# fwrite(parameter_sum_df, here(run_path, "cancer_comparison_summary.csv"))
#
