## ----setup--------------------------------------------------------------------
library("tidyr", exclude = "extract")
library("ggplot2")
library("stringr")
library("here")
library("Matrix")
library("magrittr")
library("data.table")
library("dplyr")
library("deconvR")

set.seed(123)

here::i_am("scripts/deconvolution_benchmark.R")

source(here("functions/dedupe_sigmut_mat.R"))
source(here("functions/rmse.R"))


## ----parameters --------------------------------------------------------------
testing <- TRUE

# set threshold for which transcript count is necessary for a transcript to be
# considered as an indicator for a cell type
count_thresh <- 1

n_repeat <- 10
pseudobulk_cell_frac <- 0.1

## ----data_loading-------------------------------------------------------------
data_full_meta <- fread(here(
  "datasets/Wu_etal_2021_BRCA_scRNASeq/metadata.csv"
)) %>%
  rename(cell = V1)

data_full_rownames <- readLines(here(
  "datasets/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv"
))

data_full_colnames <- readLines(here(
  "datasets/Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv"
))

data_full_matrix <- readMM(here(
  "datasets/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx"
))

data_full_matrix@Dimnames <- list(
  data_full_rownames,
  data_full_colnames
)

## ----set_data_used------------------------------------------------------------
if (testing) {
  sample_perc <- 0.001

  rnd_row_ind <- nrow(data_full_matrix) %>%
    {
      sample(x = seq_len(.), size = sample_perc * .)
    }

  rnd_col_ind <- ncol(data_full_matrix) %>%
    {
      sample(x = seq_len(.), size = sample_perc * .)
    }

  data_subsamp_matrix <- data_full_matrix[rnd_row_ind, rnd_col_ind]

  data_subsamp_meta <- data_full_meta[rnd_col_ind, ]

  count_mat <- data_subsamp_matrix
  meta <- data_subsamp_meta
} else {
  count_mat <- data_full_matrix
  meta <- data_full_meta
}

## ----convert_count_mat_to_proto_sigmat----------------------------------------
# The proto signature matrix counts how often a transcript was found in each
# cell type
celltypes <- meta %>%
  extract2("celltype_major") %>%
  unique()

celltype_cell_map <- lapply(
  celltypes,
  function(celltype, meta_df, cell_colname, celltype_colname) {
    meta_df %>%
      filter(if_any(all_of(celltype_colname), ~ .x == celltype)) %>%
      as.data.frame() %>%
      extract2(cell_colname)
  },
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
  t() %>%
  as.matrix() %>%
  rowsum(group = celltype_group) %>%
  t()


## ----signature_matrix_generation----------------------------------------------
# TODO: Consider transcript counts as weights.

sigmat <- proto_sigmat %>%
  is_weakly_greater_than(count_thresh) %>%
  # Simple as.numeric() returns a vector.
  multiply_by(1)


## ----pseudobulk_generation----------------------------------------------------
n_bulk_cells <- pseudobulk_cell_frac * ncol(count_mat)

pseudobulk_list <-
  # predraw which cells are used in each pseudobulk
  lapply(
    rep(list(seq_len(ncol(count_mat))), n_repeat),
    function(sample_vec, n) {
      sample(sample_vec, n)
    },
    n_bulk_cells
  ) %>%
  # actually draw pseudobulks
  lapply(
    function(idx_vec, count_mat, celltype_map) {
      bulk_count_mat <- count_mat[, idx_vec] %>%
        # FIXME: Decomplicate, this just removes transcripts not found in the
        # pseudo-bulk
        extract(which(rowSums(.) > 0), seq_len(ncol(.)))

      transcript_counts <- rowSums(bulk_count_mat)

      bulk_cells <- colnames(bulk_count_mat)

      celltype_counts <- celltype_map %>%
        extract(which(. %in% bulk_cells)) %>%
        names() %>%
        table()

      return(
        list(
          count_mat = bulk_count_mat,
          transcript_counts = transcript_counts,
          celltype_counts = celltype_counts
        )
      )
    },
    count_mat, celltype_cell_map
  )


## ----deconvolution------------------------------------------------------------
bulk_transcript_counts <- pseudobulk_list[[1]][["transcript_counts"]]
bulk_count_mat         <- pseudobulk_list[[1]][["count_mat"]]
bulk_celltype_counts   <- pseudobulk_list[[1]][["celltype_counts"]]
bulk_transcript_idx <- which(
  rownames(sigmat) %in% names(bulk_transcript_counts)
)

bulk_sigmat <- sigmat[bulk_transcript_idx, ]

bulk_sigmat_deduped <- dedupe_sigmut_mat(bulk_sigmat)

deconv_bulk <- bulk_transcript_counts %>%
  {
    data.frame(
      IDs = names(.),
      sample = . / sum(.)
    )
  }

deconv_reference <- bulk_sigmat_deduped %>%
  as.data.frame() %>%
  mutate(IDs = rownames(.)) %>%
  select(IDs, everything())

deconv_props <- deconvolute(
  reference = deconv_reference,
  bulk = deconv_bulk,
  model = "nnls"
)$proportions


## ----compute_deconv_err-------------------------------------------------------
true_prop_df <- bulk_celltype_counts %>%
  as.matrix() %>%
  divide_by(sum(.)) %>%
  {
    data.frame(celltype = rownames(.), prop = .)
  }

deconv_prop_df <- deconv_props %>%
  as.matrix() %>%
  extract(i = 1, j = ) %>%
  {
    data.frame(celltype = names(.), prop = .)
  }

all_prop_df <- left_join(
  true_prop_df,
  deconv_prop_df,
  by = "celltype",
  suffix = c("_true", "_deconv")
)

deconv_err <- with(all_prop_df, rmse(prop_true, prop_deconv))
deconv_err


## ----plot_deconv_err----------------------------------------------------------
all_prop_df %>%
  pivot_longer(matches("prop_*"),
    names_prefix = "prop_",
    names_to = "source",
    values_to = "prop"
  ) %>%
  ggplot(aes(celltype, prop, fill = source)) +
  geom_col(alpha = 0.5, position = position_identity()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

