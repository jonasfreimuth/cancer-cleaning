## ----setup--------------------------------------------------------------------
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

set.seed(123)

here::i_am("scripts/deconvolution_benchmark.R")

source(here("functions/dedupe_sigmut_mat.R"))
source(here("functions/rmse.R"))


## ----parameters --------------------------------------------------------------
testing <- TRUE

# Step size for exploring the effect of the count threshold, i.e. the threshold
# for which transcript count is necessary for a transcript to be considered as
# an indicator for a cell type (strict greater than, applied after
# normalization)
count_thresh_step_frac <- 0.1

n_repeat <- 200
pseudobulk_cell_frac <- 0.1

cat(paste0(
  "\nRun params:\n",
  "\tTest run: ", testing, "\n",
  "\tCount matrix threshold step size: ", count_thresh_step_frac, "\n",
  "\tNumber of repeat samplings: ", n_repeat, "\n",
  "\tFraction of ground truth sampled per pseudobulk: ", pseudobulk_cell_frac,
  "\n"
))

## ----functions ---------------------------------------------------------------
load_experiment <- function(count_mat_file, rowname_file, colname_file,
                            meta_file, testing) {
  # FIXME Make this function more general
  data_full_meta <- fread(meta_file) %>%
    rename(cell = V1)

  data_full_rownames <- readLines(rowname_file)

  data_full_colnames <- readLines(colname_file)

  # FIXME Turn this hack for saving time into something robust
  if (!exists("data_full_matrix")) {
    data_full_matrix <- readMM(count_mat_file)

    data_full_matrix@Dimnames <- list(
      data_full_rownames,
      data_full_colnames
    )
  }

  ## ----set_data_used----------------------------------------------------------
  if (testing) {
    sample_perc <- 0.01

    rnd_row_ind <- nrow(data_full_matrix) %>%
      {
        sample(x = seq_len(.), size = sample_perc * .)
      }

    rnd_col_ind <- ncol(data_full_matrix) %>%
      {
        sample(x = seq_len(.), size = sample_perc * .)
      }

    # FIXME Turn this hack for saving time into something robust
    if (!exists("data_subsamp_matrix")) {
      data_subsamp_matrix <- data_full_matrix[rnd_row_ind, rnd_col_ind]
    }

    data_subsamp_meta <- data_full_meta %>%
      filter(cell %in% colnames(data_subsamp_matrix))

    count_mat <- data_subsamp_matrix
    meta <- data_subsamp_meta
  } else {
    count_mat <- data_full_matrix
    meta <- data_full_meta
  }

  return(list(
    count_mat = count_mat,
    meta = meta
  ))
}


lognorm_row <- function(count_row, base = 10) {
  row_sum <- sum(count_row)

  # Remove rows without transcripts (might happen if down sampling)
  if (row_sum == 0) {
    return(NULL)
  }

  # log normalize counts
  count_row <- log((count_row / row_sum) + 1, base = base) * 10^4

  return(count_row)
}


seq_base <- function(start, stop, step_frac, base = 10) {
  in_vec <- c(start, stop)
  zero_ind_vec <- in_vec == 0
  if (any(zero_ind_vec)) {
    in_vec[zero_ind_vec] <- in_vec[zero_ind_vec] + 10^-16

    start <- in_vec[1]
    stop <- in_vec[2]
  }


  start_base <- log(start, base = base)
  stop_base <- log(stop, base = base)

  seq_len <- stop_base - start_base
  step_size <- seq_len * step_frac
  step_seq <- seq(from = start_base, to = stop_base, by = step_size)

  return(base^step_seq)
}


uniquify_sigmat <- function(sigmat) {
  .row_is_identifiying <- function(row, frac = 0) {
    # Determine if a row can be used to identify a column. Whether or not that
    # is the case is determined by whether a fraction of elements is non-zero.
    # The fraction maps to the number of elements which is at minimum 1 and
    # at maximum the length of the row.
    len <- length(row)
    len_out <- len * frac

    if (len_out < 1) {
      len_out <- 1
    } else if (len_out > len) {
      len_out <- len
    }

    row_sum <- sum(row > 0)

    row_sum <= len_out && row_sum > 0
  }

  unique_transcript_idx_vec <- sigmat %>%
    apply(
      1,
      .row_is_identifiying
    )

  # if only one row is selected, we would get a vector
  sigmat_unique <- sigmat[unique_transcript_idx_vec, , drop = FALSE]

  return(sigmat_unique)
}


reference_from_thresh <- function(count_thresh, proto_sigmat) {
  sigmat <- proto_sigmat %>%
    is_greater_than(count_thresh) %>%
    # Simple as.numeric() returns a vector.
    multiply_by(1) %>%
    uniquify_sigmat() %>%
    dedupe_sigmut_mat()

  # sigmats with colSums equal 0 lead to deconv troubles
  # TODO check what problems this causes, e.g. leading to
  # 1 or 0 col sigmats.
  col_sums <- colSums(sigmat)

  sigmat <- sigmat[, col_sums > 0, drop = FALSE]

  if (!ncol(sigmat) < nrow(sigmat)) {
    warning(paste0(
      "Threshold ", count_thresh, " produced a non-tall signature matrix, ",
      "NULL returned instead. (Will likely be filtered out downstream.)"
    ))
    return(NULL)
  }

  deconv_ref <- sigmat %>%
    as.data.frame() %>%
    mutate(IDs = rownames(.)) %>%
    select(IDs, everything())

  return(deconv_ref)
}


create_celltype_map <- function(celltype, meta_df, cell_colname,
                                celltype_colname) {
  meta_df %>%
    filter(if_any(all_of(celltype_colname), ~ .x == celltype)) %>%
    as.data.frame() %>%
    extract2(cell_colname)
}


pseudobulk_from_idx <- function(idx_vec, count_mat, celltype_map) {
  bulk_count_mat <- count_mat[, idx_vec] %>%
    scuttle::calculateTPM()

  transcript_counts <- rowSums(bulk_count_mat)

  bulk_cells <- colnames(bulk_count_mat)

  celltype_counts <- celltype_map %>%
    extract(which(. %in% bulk_cells)) %>%
    names() %>%
    table()

  return(
    list(
      transcript_counts = transcript_counts,
      celltype_counts = celltype_counts
    )
  )
}


deconvolute_pseudobulk <- function(pseudobulk, deconv_ref) {
  transcript_props <- pseudobulk[["transcript_counts"]] %>%
    extract(names(.) %in% deconv_ref$IDs) %>%
    divide_by(sum(.))
  celltype_props <- pseudobulk[["celltype_counts"]] %>%
    as.matrix() %>%
    divide_by(sum(.))

  deconv_bulk <- transcript_props %>%
    {
      data.frame(
        IDs = names(.),
        sample = .
      )
    }

  capture.output(
    suppressMessages(
      deconv_props <- deconvolute(
        reference = deconv_ref,
        bulk = deconv_bulk,
        model = "qp"
      )$proportions
    ),
    type = c("output")
  )

  sigmat <- deconv_ref %>%
    select(-IDs) %>%
    as.matrix()

  dp_vec <- deconv_props %>%
    as.matrix() %>%
    extract(i = 1, j = )

  transcript_props_pred <- sigmat %*% dp_vec

  deconv_resid <- transcript_props - transcript_props_pred

  true_prop_df <- celltype_props %>%
    {
      data.frame(celltype = rownames(.), prop = .)
    } %>%
    separate_rows(celltype, sep = "_")

  deconv_prop_df <- dp_vec %>%
    {
      data.frame(celltype = names(.), prop = .)
    }

  deconv_res <- full_join(
    true_prop_df,
    deconv_prop_df,
    by = "celltype",
    suffix = c("_true", "_deconv")
  ) %>%
    # ctypes not found in deconv output (likely due to them not having cols
    # uniquely identifying them in the deconv_ref) should be set to 0 to be able
    # to still compute errors properly.
    # TODO Does this make sense? Is this comparable to using Others col, and
    # should that maybe be preferred visavis cancer prop calculation from
    # residuals?
    mutate(across(matches("prop_*"), ~ replace_na(.x, 0))) %>%
    mutate(
      abs_err = abs(prop_true - prop_deconv),
      n_ctypes = nrow(deconv_prop_df)
    )

  return(list(
    res = deconv_res,
    residuals = deconv_resid
  ))
}


benchmark_reference <- function(deconv_ref, pseudobulk_list) {
  deconv_res_list <- pseudobulk_list %>%
    lapply(
      deconvolute_pseudobulk,
      deconv_ref
    )

  deconv_prop_list <- deconv_res_list %>%
    lapply(function(deconv_res) {
      deconv_res$res
    })

  deconv_err_vec <- deconv_prop_list %>%
    lapply(
      function(deconv_prop_df) {
        with(deconv_prop_df, rmse(prop_true, prop_deconv))
      }
    ) %>%
    unlist()

  print(mean(deconv_err_vec))

  all_prop_df <- deconv_prop_list %>%
    bind_rows(.id = "sample")

  all_prop_sum_df <- all_prop_df %>%
    group_by(sample) %>%
    summarise(
      across(all_of(c(
        "prop_true", "prop_deconv"
      )), sum),
      med_abs_err = median(abs_err),
      unq_n_ct = paste(unique(n_ctypes), sep = "_"),
      rmse = rmse(prop_true, prop_deconv)
    )

  return(list(
    "errors" = deconv_err_vec,
    "deconv_res" = all_prop_df,
    "deconv_sum" = all_prop_sum_df
  ))
}


plot_deconv_res <- function(deconv_res) {
  all_prop_df <- deconv_res[["deconv_res"]]
  deconv_err_vec <- deconv_res[["errors"]]


  all_prop_df %>%
    group_by(celltype) %>%
    summarise(
      rmse = rmse(prop_true, prop_deconv)
    ) %>%
    ggplot(aes(celltype, rmse)) +
    geom_col(alpha = 0.5, position = position_identity()) +
    labs(
      title = paste(
        "Mean sample RMSE:",
        mean(deconv_err_vec) %>% round(2)
      ),
      x = "Cell types",
      y = "Per celltype RMSE"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
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

count_mat <- data$count_mat
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
  scuttle::calculateTPM() %>%
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
#   set_names(as.character(.))

# temp solution
count_thresh_vec <- 7323

## ----signature_matrix_generation----------------------------------------------
# TODO: Consider transcript counts as weights.
# TODO: Explore effects of sigmat threshold
# TODO: Add others?
# TODO: Do the residual computation and correlate with actual cancer prop

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
    count_mat, celltype_cell_map
  )


## ----deconvolution------------------------------------------------------------
deconv_res_list <- lapply(
  deconv_ref_list,
  benchmark_reference,
  pseudobulk_list
)


## ----plot_deconv_err----------------------------------------------------------
plot_list <- lapply(
  deconv_res_list,
  plot_deconv_res
)
