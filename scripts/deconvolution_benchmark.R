## ----setup--------------------------------------------------------------------
library("cowplot")
library("scuttle")
library("tidyr", exclude = "extract")
library("ggplot2")
library("stringr")
library("here")
library("Matrix")
library("magrittr")
library("data.table")
library("dplyr")
library("deconvR")

here::i_am("scripts/deconvolution_benchmark.R")

source(here("functions/dedupe_sigmut_mat.R"))
source(here("functions/rmse.R"))

source(here("functions/util_functions.R"))
source(here("functions/norm_functions.R"))
source(here("functions/benchmark_functions.R"))

source(here("etc/plot_themes.R"))

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

  # Whether the signature matrix is supposed to be binary. This is achieved
  # via thresholding along a equally spaced exponential sequence.
  # There are also some finer difference how feature selection is performed.
  "binary_sigmat",

  # Step size for exploring the effect of the count threshold, i.e. the
  # threshold for which transcript count is necessary for a transcript to be
  # considered as an indicator for a cell type (strict greater than, applied
  # after normalization)
  # Has no effect if a non-binary signature matrix is requested.
  "count_thresh_step_frac",

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
  # * (tpm & logtpm: Deprecated.)
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
  binary_sigmat = "FALSE",
  count_thresh_step_frac = "0.3",
  n_pseudobulk = "10",
  pseudobulk_cell_frac = "0.2",
  normalization_type = "lognorm",
  normalize_independently = "TRUE",
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
params$normalize_independently %<>%
  as.logical()
params$binary_sigmat %<>%
  as.logical()
params$count_thresh_step_frac %<>%
  as.numeric()
params$n_pseudobulk %<>%
  as.numeric()
params$pseudobulk_cell_frac %<>%
  as.numeric()
params$seed %<>%
  as.numeric()

params$norm_scale <- 10^6

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
  paste0("binary_sigmat", pair_sep, params$binary_sigmat),
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
  paste(collapse = "\n")

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


## ----data_loading-------------------------------------------------------------
data <- load_experiment(
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

count_mat <- data$count_mat %>%
  normalize_count_mat(
    type = super_norm_type,
    scale = params$norm_scale
  )
meta <- data$meta

## ----convert_count_mat_to_proto_sigmat----------------------------------------
# The proto signature matrix counts how often a transcript was found in each
# cell type
celltypes <- meta %>%
  extract2("celltype_major") %>%
  unique()

celltype_cell_map <- lapply(
  celltypes,
  create_celltype_map,
  meta, "cell", "celltype_major"
) %>%
  set_names(celltypes) %>%
  unlist() %>%
  set_names(str_extract(names(.), "[^0-9]+")) %>%
  factor(levels = colnames(count_mat)) %>%
  sort()

celltype_group <- celltype_cell_map %>%
  names()

proto_sigmat <- count_mat %>%
  t() %>%
  as.matrix() %>%
  # TODO Should cols instead be averaged?
  # The case against would be that in bulks the number of celltypes has an
  # effect on the expression values (i.e. more cells in bulk -> higher
  # transcript counts)
  rowsum(group = celltype_group) %>%
  t() %>%
  normalize_count_mat(
    type = sub_norm_type,
    scale = params$norm_scale
  )

count_mat <- count_mat %>%
  as.matrix()


## ----signature_matrix_generation----------------------------------------------
if (params$binary_sigmat) {
  count_range <- proto_sigmat %>%
    as.vector() %>%
    extract(. > 0) %>%
    range()

  count_thresh_vec <- seq_base(
    count_range[1],
    count_range[2],
    params$count_thresh_step_frac,
    base = 3
  ) %>%
    {
      c(0, .)
    } %>%
    set_names(as.character(round(., 2)))

  deconv_ref_list <- lapply(
    count_thresh_vec,
    bin_reference_from_thresh,
    proto_sigmat = proto_sigmat
  )
} else {
  count_thresh_vec <- 0 %>%
    set_names(as.character(round(., 2)))

  deconv_ref_list <- proto_sigmat %>%
    reference_non_binary() %>%
    list() %>%
    set_names("0")
}

# For the moment, the reference needs to always be a list of length 1 with the
# name of the element
is_null_reference <- lapply(deconv_ref_list, is.null) %>%
  unlist()

deconv_ref_list <- deconv_ref_list %>%
  extract(!is_null_reference)

## ----pseudobulk_generation----------------------------------------------------
n_bulk_cells <- params$pseudobulk_cell_frac * ncol(count_mat)

pseudobulk_list <-
  # predraw which cells are used in each pseudobulk
  lapply(
    rep(list(seq_len(ncol(count_mat))), params$n_pseudobulk),
    sample,
    n_bulk_cells
  ) %>%
  # actually draw pseudobulks
  lapply(
    pseudobulk_from_idx,
    count_mat, celltype_cell_map,
    norm_type = sub_norm_type,
    scale = params$norm_scale
  )


## ----deconvolution------------------------------------------------------------
split_vals <- c(TRUE, FALSE)

split_res_list <- lapply(
  split_vals,
  function(split) {
    lapply(
      deconv_ref_list,
      benchmark_reference,
      pseudobulk_list,
      split_cancer = split,
      deconv_method = params$deconv_method
    )
  }
) %>%
  set_names(split_vals)

## ----plot_deconv_err----------------------------------------------------------
sample_sum_df <- split_res_list %>%
  dfextract2("deconv_sum", "sigmat_thresh", "split")

rename_vec <- sample_sum_df %>%
  head(n = 1) %>%
  select(matches(c("rmse", "cancer_expr_v_*"))) %>%
  names()

rename_targets <- rename_vec %>%
  str_replace_all(pattern = "cancer_expr", replacement = "cexpr") %>%
  str_replace_all(pattern = "rmse", replacement = "bulk_rmse")

rename_vec <- rename_vec %>%
  set_names(rename_targets)

parameter_sum_df <- sample_sum_df %>%
  group_by(split, sigmat_thresh) %>%
  rename(!!!rename_vec) %>%
  summarize(
    across(
      matches(c("bulk_rmse", "cexpr_v_*")),
      list("mean" = mean),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop_last"
  )


celltype_rmse_df <- split_res_list %>%
  dfextract2("deconv_res", "sigmat_thresh", "split")

celltype_df_sum <- celltype_rmse_df %>%
  group_by(split, sigmat_thresh, celltype) %>%
  summarise(
    # per celltype rmse, across all samples
    celltype_rmse = rmse(prop_true, prop_deconv),
    .groups = "drop_last"
  ) %>%
  left_join(parameter_sum_df, by = c("split", "sigmat_thresh")) %>%
  mutate(sigmat_thresh = as.numeric(sigmat_thresh))


plot_err <- ggplot(
  celltype_df_sum,
  aes(celltype, celltype_rmse)
) +
  geom_col(alpha = 0.5, position = position_identity()) +
  geom_text(
    aes(
      x = length(unique(celltype)) / 2,
      y = max(celltype_rmse) * 1.1,
      label = paste0(
        "Mean bulk RMSE: ", round(mean_bulk_rmse, 3)
      )
    ),
    vjust = 1
  ) +
  labs(
    x = "Cell types",
    y = "Per celltype RMSE"
  ) +
  facet_grid(
    cols = vars(split),
    rows = vars(sigmat_thresh)
  ) +
  theme_benchmark() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dir.create(here(run_path, "plots"), recursive = TRUE, showWarnings = FALSE)

ggsave(here(run_path, "plots", "rmse_plot.png"), plot_err,
  height = params$base_height + (
    length(count_thresh_vec) * params$facet_height
  ),
  width = params$base_width + (
    length(split_vals) * params$facet_width
  )
)

# Resid vs Expr Plot ------------------------------------------------------
resid_expr_marker_df <- split_res_list %>%
  dfextract2("resid_expr_marker_df", "sigmat_thresh", "split") %>%
  left_join(parameter_sum_df, by = c("split", "sigmat_thresh"))

plot_resid_expr_marker <- ggplot(
  resid_expr_marker_df,
  aes(cancer_expr, resid, col = sample)
) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_text(
    aes(
      x = mean(c(min(cancer_expr), max(cancer_expr))),
      y = max(resid) * 1.1,
      label = paste0(
        "Mean r: ", round(mean_cexpr_v_resid_marker, 3)
      )
    ),
    vjust = 1,
    col = "gray33"
  ) +
  labs(
    x = "Cancer expression",
    y = "Transcript residual"
  ) +
  facet_wrap(~sample) +
  facet_grid(
    cols = vars(split),
    rows = vars(sigmat_thresh),
    scales = "free"
  ) +
  theme_benchmark() +
  theme(legend.position = "none")

dir.create(here(run_path, "plots"), recursive = TRUE, showWarnings = FALSE)

ggsave(here(run_path, "plots", "resid_expr_marker_plot.png"),
  plot_resid_expr_marker,
  height = params$base_height + (
    length(count_thresh_vec) * params$facet_height
  ),
  width = params$base_width + (
    length(split_vals) * params$facet_width
  )
)

resid_expr_all_df <- split_res_list %>%
  dfextract2("resid_expr_all_df", "sigmat_thresh", "split") %>%
  left_join(parameter_sum_df, by = c("split", "sigmat_thresh"))

plot_resid_expr_all <- ggplot(
  resid_expr_all_df,
  aes(cancer_expr, resid, col = sample)
) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_text(
    aes(
      x = mean(c(min(cancer_expr), max(cancer_expr))),
      y = max(resid) * 1.1,
      label = paste0(
        "Mean r: ", round(mean_cexpr_v_resid_all, 3)
      )
    ),
    vjust = 1,
    col = "gray33"
  ) +
  labs(
    x = "Cancer expression",
    y = "Transcript residual"
  ) +
  facet_wrap(~sample) +
  facet_grid(
    cols = vars(split),
    rows = vars(sigmat_thresh),
    scales = "free"
  ) +
  theme_benchmark() +
  theme(legend.position = "none")

dir.create(here(run_path, "plots"), recursive = TRUE, showWarnings = FALSE)

ggsave(here(run_path, "plots", "resid_expr_all_plot.png"), plot_resid_expr_all,
  height = params$base_height + (
    length(count_thresh_vec) * params$facet_height
  ),
  width = params$base_width + (
    length(split_vals) * params$facet_width
  )
)


# Pred Trans Abund vs Expr Plot -------------------------------------------
plot_pred_prop_expr_marker <- ggplot(
  resid_expr_marker_df,
  aes(cancer_expr, deconv_pred, col = sample)
) +
  geom_point() +
  geom_text(
    aes(
      x = mean(c(min(cancer_expr), max(cancer_expr))),
      y = max(deconv_pred) * 1.2,
      label = paste0(
        "Mean r: ", round(mean_cexpr_v_deconv_pred_marker, 3)
      )
    ),
    vjust = 1,
    col = "gray33"
  ) +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    x = "Cancer expression",
    y = "Deconvolution model prediction"
  ) +
  facet_wrap(~sample) +
  facet_grid(
    cols = vars(split),
    rows = vars(sigmat_thresh),
    scales = "free"
  ) +
  theme_benchmark() +
  theme(legend.position = "none")

ggsave(here(run_path, "plots", "pred_prop_expr_marker_plot.png"),
  plot_pred_prop_expr_marker,
  height = params$base_height + (
    length(count_thresh_vec) * params$facet_height
  ),
  width = params$base_width + (
    length(split_vals) * params$facet_width
  )
)


plot_pred_prop_expr_all <- ggplot(
  resid_expr_all_df,
  aes(cancer_expr, deconv_pred, col = sample)
) +
  geom_point() +
  geom_text(
    aes(
      x = mean(c(min(cancer_expr), max(cancer_expr))),
      y = max(deconv_pred) * 1.2,
      label = paste0(
        "Mean r: ", round(mean_cexpr_v_deconv_pred_all, 3)
      )
    ),
    vjust = 1,
    col = "gray33"
  ) +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    x = "Cancer expression",
    y = "Deconvolution model prediction"
  ) +
  facet_wrap(~sample) +
  facet_grid(
    cols = vars(split),
    rows = vars(sigmat_thresh),
    scales = "free"
  ) +
  theme_benchmark() +
  theme(legend.position = "none")

ggsave(here(run_path, "plots", "pred_prop_expr_all_plot.png"),
  plot_pred_prop_expr_all,
  height = params$base_height + (
    length(count_thresh_vec) * params$facet_height
  ),
  width = params$base_width + (
    length(split_vals) * params$facet_width
  )
)

# ----Data saving---------------------------------------------------------------
file.copy(
  here("scripts/deconvolution_benchmark.R"),
  here(run_path, "deconvolution_benchmark.R")
) %>%
  invisible()

file.copy(
  here("functions/benchmark_functions.R"),
  here(run_path, "benchmark_functions.R")
) %>%
  invisible()

print(parameter_sum_df)

fwrite(parameter_sum_df, here(run_path, "cancer_comparison_summary.csv"))
