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

source(here("functions/benchmark_functions.R"))

source(here("etc/plot_themes.R"))

## ----parameters --------------------------------------------------------------
if (!exists("script_args")) {
  script_args <- commandArgs(trailingOnly = TRUE)
}

# names must match script params
arg_names <- c(
  "data_path",
  "count_thresh_step_frac",
  "n_repeat",
  "pseudobulk_cell_frac",
  "normalization_type",
  "normalize_independently",
  "deconv_method",
  "seed"
)

if (length(script_args) > 0) {
  cat("\nUsing provided args...\n")
  names(script_args) <- arg_names
} else {
  cat("\nUsing default args...\n")
  script_args <- c(
    data_path = "datasets/Wu_etal_downsampled_test/",

    # Step size for exploring the effect of the count threshold, i.e. the
    # threshold for which transcript count is necessary for a transcript to be
    # considered as an indicator for a cell type (strict greater than, applied
    # after normalization)
    count_thresh_step_frac = "0.3",
    n_repeat = "200",
    pseudobulk_cell_frac = "0.2",
    normalization_type = "tpm",
    normalize_independently = "TRUE",
    deconv_method = "qp",
    seed = "123"
  )
}

for (i in seq_along(script_args)) {
  # FIXME Change the script to use args vector and get rid of this assign hack.
  assign(names(script_args[i]), script_args[i])
}

normalize_independently %<>%
  as.logical()
count_thresh_step_frac %<>%
  as.numeric()
n_repeat %<>%
  as.numeric()
pseudobulk_cell_frac %<>%
  as.numeric()
seed %<>%
  as.numeric()


base_width <- 3
base_height <- 2

facet_height <- 2
facet_width <- 5.5

param_sep <- "_"
pair_sep <- "-"

parameter_string <- paste(
  paste0("data_path", pair_sep, basename(data_path)),
  paste0("seed", pair_sep, seed),
  paste0("method", pair_sep, deconv_method),
  paste0("indepnorm", pair_sep, normalize_independently),
  paste0("normtype", pair_sep, normalization_type),
  paste0("nrepeat", pair_sep, n_repeat),
  paste0("sizefrac", pair_sep, pseudobulk_cell_frac),
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

run_info <- paste0(
  "Run params:\n",
  "\tDataset used:\t", data_path, "\n",
  "\tRandomization seed:\t", seed, "\n",
  "\tDeconvolution method:\t", deconv_method, "\n",
  "\tIndependent normalization of count matrix for bulk and reference:\t",
  normalize_independently, "\n",
  "\tCount matrix threshold step size:\t", count_thresh_step_frac, "\n",
  "\tNumber of repeat samplings:\t", n_repeat, "\n",
  "\tFraction of ground truth sampled per pseudobulk:\t", pseudobulk_cell_frac,
  "\n", "\n",
  "Outputs will be found at ", run_path, "\n"
)

cat(run_info, file = here(run_path, "run_info.txt"))

cat(run_info)

set.seed(seed)

if (normalize_independently) {
  super_norm_type <- NULL
  sub_norm_type <- normalization_type
} else {
  super_norm_type <- normalization_type
  sub_norm_type <- NULL
}


## ----data_loading-------------------------------------------------------------
data <- load_experiment(
  count_mat_file = here(
    data_path, "count_matrix_sparse.mtx"
  ),
  rowname_file = here(
    data_path, "count_matrix_genes.tsv"
  ),
  colname_file = here(
    data_path, "count_matrix_barcodes.tsv"
  ),
  meta_file = here(
    data_path, "metadata.csv"
  )
)

count_mat <- data$count_mat %>%
  normalize_count_mat(type = super_norm_type)
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
  normalize_count_mat(type = sub_norm_type) %>%
  t() %>%
  as.matrix() %>%
  rowsum(group = celltype_group) %>%
  t()

count_mat <- count_mat %>%
  as.matrix()

## ----signature_matrix_generation----------------------------------------------
# TODO: Consider transcript counts as weights.
# TODO: Add Others col?
count_range <- proto_sigmat %>%
  as.vector() %>%
  extract(. > 0) %>%
  range()

count_thresh_vec <- seq_base(
  count_range[1],
  count_range[2],
  count_thresh_step_frac,
  base = 3
) %>%
  {
    c(0, .)
  } %>%
  set_names(as.character(round(., 2)))

deconv_ref_list <- lapply(
  count_thresh_vec,
  reference_from_thresh,
  proto_sigmat = proto_sigmat
)

is_null_reference <- lapply(deconv_ref_list, is.null) %>%
  unlist()

deconv_ref_list <- deconv_ref_list %>%
  extract(!is_null_reference)

## ----pseudobulk_generation----------------------------------------------------
n_bulk_cells <- pseudobulk_cell_frac * ncol(count_mat)

pseudobulk_list <-
  # predraw which cells are used in each pseudobulk
  lapply(
    rep(list(seq_len(ncol(count_mat))), n_repeat),
    sample,
    n_bulk_cells
  ) %>%
  # actually draw pseudobulks
  lapply(
    pseudobulk_from_idx,
    count_mat, celltype_cell_map,
    norm_type = sub_norm_type
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
      deconv_method = deconv_method
    )
  }
) %>%
  set_names(split_vals)

## ----plot_deconv_err----------------------------------------------------------
sample_sum_df <- split_res_list %>%
  dfextract2("deconv_sum", "sigmat_thresh", "split")

parameter_sum_df <- sample_sum_df %>%
  group_by(split, sigmat_thresh) %>%
  summarize(
    mean_bulk_rmse = mean(rmse),
    mean_canc_expr_corr = mean(cancer_expr_corr),
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
  geom_text(aes(
    x = length(unique(celltype)) / 2,
    y = max(celltype_rmse) * 1.1,
    label = paste0(
      "Mean bulk RMSE: ", round(mean_bulk_rmse, 3), ",\n",
      "Mean r: ", round(mean_canc_expr_corr, 3)
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
  height = base_height + (length(count_thresh_vec) * facet_height),
  width = base_width + (length(split_vals) * facet_width)
)

# Resid vs Expr Plot ------------------------------------------------------
resid_expr_df <- split_res_list %>%
  dfextract2("resid_expr_df", "sigmat_thresh", "split")

plot_resid_expr <- ggplot(
  resid_expr_df,
  aes(prop, resid, col = sample)
) +
  geom_point() +
  labs(
    x = "Transcript Proportion",
    y = "Transcript Residual"
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

ggsave(here(run_path, "plots", "resid_expr_plot.png"), plot_resid_expr,
  height = base_height + (length(count_thresh_vec) * facet_height),
  width = base_width + (length(split_vals) * facet_width)
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
