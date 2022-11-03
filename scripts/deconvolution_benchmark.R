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

## ----parameters --------------------------------------------------------------
testing <- TRUE

# Step size for exploring the effect of the count threshold, i.e. the threshold
# for which transcript count is necessary for a transcript to be considered as
# an indicator for a cell type (strict greater than, applied after
# normalization)
count_thresh_step_frac <- 0.1

n_repeat <- 200
pseudobulk_cell_frac <- 0.1
normalization_type <- "tpm"
normalize_independently <- TRUE
deconv_method <- "qp"

seed <- 123

base_width <- 3
base_height <- 2

facet_height <- 1
facet_width <- 4.5

param_sep <- "_"
pair_sep <- "-"

parameter_string <- paste(
  paste0("testing", pair_sep, testing),
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
    format(Sys.time(), "%Y%m%d-%H%M"),
    parameter_string,
    sep = param_sep
  )
)

dir.create(run_path, recursive = TRUE, showWarnings = FALSE)

cat(paste0(
  "\nRun params:\n",
  "\tTest run: ", testing, "\n",
  "\tCount matrix threshold step size: ", count_thresh_step_frac, "\n",
  "\tNumber of repeat samplings: ", n_repeat, "\n",
  "\tFraction of ground truth sampled per pseudobulk: ", pseudobulk_cell_frac,
  "\n", "\n",
  "Outputs will be found at ", run_path, "\n"
))

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
    "datasets/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx"
  ),
  rowname_file = here(
    "datasets/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv"
  ),
  colname_file = here(
    "datasets/Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv"
  ),
  meta_file = here(
    "datasets/Wu_etal_2021_BRCA_scRNASeq/metadata.csv"
  ),
  testing = testing
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
  set_names(str_extract(names(.), "[^0-9]+"))

celltype_group <- celltype_cell_map %>%
  factor(levels = colnames(count_mat)) %>%
  sort() %>%
  names()

proto_sigmat <- count_mat %>%
  normalize_count_mat(type = sub_norm_type) %>%
  t() %>%
  as.matrix() %>%
  rowsum(group = celltype_group) %>%
  t()

count_mat <- count_mat %>%
  as.matrix()

# count_range <- proto_sigmat %>%
#   as.vector() %>%
#   extract(. > 0) %>%
#   range()
#
# count_thresh_vec <- seq_base(
#   count_range[1],
#   count_range[2],
#   count_thresh_step_frac,
#   base = 3
# ) %>%
#   {
#     c(0, .)
#   } %>%
#   set_names(as.character(round(., 2)))

# temp solution
count_thresh_vec <- c(7323) %>%
  set_names(as.character(.))

## ----signature_matrix_generation----------------------------------------------
# TODO: Consider transcript counts as weights.
# TODO: Add Others col?

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
    count_mat, celltype_cell_map, sub_norm_type
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

mean_rmses <- split_res_list %>%
  lapply(
    function(split_res) {
      split_res %>%
        lapply(
          function(deconv_res) {
            deconv_res %>%
              extract2("errors") %>%
              mean() %>%
              {
                data.frame(mean_rmse = .)
              }
          }
        ) %>%
        bind_rows(.id = "sigmat_thresh")
    }
  ) %>%
  bind_rows(.id = "split")

cancer_comp_df <- split_res_list %>%
  lapply(
    function(split_res) {
      split_res %>%
        lapply(
          function(deconv_res) {
            deconv_res %>%
              extract2("cancer_comp")
          }
        ) %>%
        bind_rows(.id = "sigmat_thresh")
    }
  ) %>%
  bind_rows(.id = "split")

cor_meth <- "pearson"

cancer_comp_df_sum <- cancer_comp_df %>%
  group_by(split, sigmat_thresh) %>%
  summarize(
    cor_resid_by_sigmat = cor(prop_true, by_sigmat, method = cor_meth),
    cor_resid_by_transcr_prop = cor(
      prop_true, by_transcr_prop,
      method = cor_meth
    ),
    cor_sum_sq_res = cor(prop_true, sum_sq_resid, method = cor_meth),
    cor_sum_abs_res = cor(prop_true, sum_abs_resid, method = cor_meth),
    cor_sum_res = cor(prop_true, sum_resid, method = cor_meth)
  ) %>%
  left_join(mean_rmses, by = c("split", "sigmat_thresh"))

print(cancer_comp_df_sum)

corr_cols <- c(
  "by_sigmat" = "cor_resid_by_sigmat",
  "by_transcr_prop" = "cor_resid_by_transcr_prop",
  "sum_sq_resid" = "cor_sum_sq_res",
  "sum_abs_resid" = "cor_sum_abs_res",
  "sum_resid" = "cor_sum_res"
)

corr_col_info <- data.frame(
  col = names(corr_cols),
  corr = corr_cols,
  display_name = c(
    "Residuals *\ncancer signature",
    "Residuals *\ncancer transcript proportions",
    "Sum of squared\nresiduals",
    "Sum of absolute\nresiduals",
    "Sum of\nresiduals"
  )
)


plot_corr <- apply(
  corr_col_info,
  1,
  function(col_info, df, df_sum, meth) {
    rename_vec <- c("corr_col" = col_info[["col"]])
    corr_df <- df %>%
      select(
        split, sigmat_thresh, prop_true, all_of(rename_vec)
      )

    prefix <- ""

    if (cor_meth %in% c("spearman", "kendall")) {
      corr_df <- corr_df %>%
        mutate(across(
          all_of(c("prop_true", "corr_col")),
          ~ rank(.x)
        ))

      prefix <- "Rank of "
    }

    rename_vec_sum <- c("corr_col_sum" = col_info[["corr"]])
    corr_df_sum <- df_sum %>%
      select(
        split, sigmat_thresh, all_of(rename_vec_sum)
      )

    corr_df %>%
      left_join(corr_df_sum, by = c("sigmat_thresh", "split")) %>%
      mutate(sigmat_thresh = as.numeric(sigmat_thresh)) %>%
      ggplot() +
      geom_point(aes(prop_true, corr_col)) +
      geom_text(
        aes(
          x = mean(range(prop_true)),
          y = max(corr_col) + diff(range(corr_col) * 0.1),
          label = round(corr_col_sum, 3)
        )
      ) +
      labs(
        x = paste(prefix, "True cancer proportion"),
        y = paste(prefix, col_info[["display_name"]])
      ) +
      facet_grid(
        cols = vars(split),
        rows = vars(sigmat_thresh)
      ) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line()
      )
  },
  cancer_comp_df,
  cancer_comp_df_sum,
  corr_meth
) %>%
  {
    cowplot::plot_grid(
      plotlist = .,
      nrow = length(.)
    )
  }

dir.create(here(run_path, "plots"), recursive = TRUE, showWarnings = FALSE)

ggsave(here(run_path, "plots", "benchmark_plot.png"), plot_corr,
  height = base_height +
    (length(count_thresh_vec) * facet_height * nrow(corr_col_info)),
  width = base_width + (length(split_vals) * facet_width)
)

## ----plot_deconv_err----------------------------------------------------------
celltype_rmse_df <- split_res_list %>%
  lapply(
    function(split_res) {
      split_res %>%
        lapply(
          function(deconv_res) {
            mean_rmse <- deconv_res %>%
              extract2("errors") %>%
              mean()

            deconv_res %>%
              extract2("deconv_res") %>%
              group_by(celltype) %>%
              summarise(
                rmse = rmse(prop_true, prop_deconv),
                .groups = "drop_last"
              ) %>%
              mutate(overall_rmse = mean_rmse)
          }
        ) %>%
        bind_rows(.id = "sigmat_thresh")
    }
  ) %>%
  bind_rows(.id = "split")

plot_err <- celltype_rmse_df %>%
  mutate(sigmat_thresh = as.numeric(sigmat_thresh)) %>%
  ggplot(aes(celltype, rmse)) +
  geom_col(alpha = 0.5, position = position_identity()) +
  geom_text(aes(
    x = length(unique(celltype)) / 2,
    y = max(rmse) * 1.1,
    label = round(overall_rmse, 3)
  )) +
  labs(
    x = "Cell types",
    y = "Per celltype RMSE"
  ) +
  facet_grid(
    cols = vars(split),
    rows = vars(sigmat_thresh)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line()
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dir.create(here(run_path, "plots"), recursive = TRUE, showWarnings = FALSE)

ggsave(here(run_path, "plots", "rmse_plot.png"), plot_err,
  height = base_height + (length(count_thresh_vec) * facet_height),
  width = base_width + (length(split_vals) * facet_width)
)
