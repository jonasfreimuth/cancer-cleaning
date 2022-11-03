library("BiocGenerics")
library("base")
library("dplyr")
library("S4Vectors")
library("stats")
library("tidyr")
library("IRanges")
library("Matrix")
library("SummarizedExperiment")
library("GenomicRanges")
library("GenomeInfoDb")
library("deconvR")
library("magrittr")
library("scuttle")
library("utils")

# TODO Once it is clear which functions are actually useful, add proper
# docstrings.

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


normalize_count_mat <- function(count_mat, type = "lognorm", ...) {
  if (type == "lognorm") {
    norm_mat <- count_mat %>%
      apply(
        2,
        lognorm_row,
        simplify = TRUE,
        ...
      )
  } else if (type == "tpm") {
    norm_mat <- count_mat %>%
      scuttle::calculateTPM(...)
  } else {
    stop(paste0(
      "Normalization type \"", type, "\" is not supported."
    ))
  }

  return(count_mat)
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


pseudobulk_from_idx <- function(idx_vec, count_mat, celltype_map,
                                norm_type = "tpm") {
  bulk_count_mat <- count_mat[, idx_vec] %>%
    normalize_count_mat(type = norm_type)

  transcript_counts <- rowSums(bulk_count_mat)

  bulk_cells <- colnames(bulk_count_mat)

  celltype_counts <- celltype_map %>%
    extract(which(. %in% bulk_cells)) %>%
    names() %>%
    table()

  bulk_celltypes <- names(celltype_counts)

  # basically just conversion from table to named vector
  celltype_counts <- celltype_counts %>%
    as.vector() %>%
    set_names(bulk_celltypes)

  return(
    list(
      transcript_counts = transcript_counts,
      celltype_counts = celltype_counts
    )
  )
}


deconvolute_pseudobulk <- function(pseudobulk, deconv_ref,
                                   split_cancer = FALSE,
                                   cancer_cols = "Cancer Epithelial",
                                   method = "qp") {
  # FIXME Think about what multiple cancer cols would mean.
  if (length(cancer_cols) > 1) {
    stop(paste(
      "Currently more than one column of cancer celltypes (cancer_col) is",
      "not supported."
    ))
  }

  transcript_props <- pseudobulk[["transcript_counts"]] %>%
    extract(names(.) %in% deconv_ref$IDs) %>%
    divide_by(sum(.)) %>%
    extract(!is.na(.))
  celltype_props <- pseudobulk[["celltype_counts"]] %>%
    divide_by(sum(.)) %>%
    extract(!is.na(.))

  input_list <- list(
    tp = transcript_props,
    cp = celltype_props
  )

  input_lens <- lapply(input_list, length)

  if (any(input_lens == 0)) {
    return(NULL)
  }

  deconv_bulk <- transcript_props %>%
    {
      data.frame(
        IDs = names(.),
        sample = .
      )
    }

  cancer_ref <- deconv_ref %>%
    select(IDs, all_of(cancer_cols))

  if (split_cancer) {
    deconv_ref <- deconv_ref %>%
      select(-all_of(cancer_cols))
  }

  if (ncol(deconv_ref) < 3) {
    # Some deconvolution methods will fail when only a single celltype is
    # present. This ensures the correct output.
    deconv_props <- 1 %>%
      set_names(names(
        deconv_ref %>%
          select(-IDs)
      ))
  } else {
    capture.output(
      suppressMessages(
        deconv_props <- deconvolute(
          reference = deconv_ref,
          bulk = deconv_bulk,
          model = method
        )$proportions
      ),
      type = c("output")
    )
  }

  sigmat <- deconv_ref %>%
    select(-IDs) %>%
    as.matrix()

  dp_vec <- deconv_props %>%
    as.matrix() %>%
    extract(i = 1, j = )

  transcript_props_pred <- sigmat %*% dp_vec

  deconv_resid <- transcript_props - transcript_props_pred

  cancer_sig <- cancer_ref %>%
    select(all_of(cancer_cols)) %>%
    as.matrix()

  resid_cancer_prop <- t(deconv_resid) %*% cancer_sig

  cancer_trancsr_prop <- proto_sigmat %>%
    extract(
      i = rownames(.) %in% rownames(deconv_resid),
      j = colnames(.) %in% cancer_cols
    ) %>%
    divide_by(sum(.))
  prop_cancer_prop <- t(deconv_resid) %*% cancer_trancsr_prop

  cancer_prop <- data.frame(
    by_sigmat = as.vector(resid_cancer_prop),
    by_transcr_prop = as.vector(prop_cancer_prop)
  )

  true_prop_df <- celltype_props %>%
    {
      data.frame(celltype = names(.), prop = .)
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
    residuals = deconv_resid,
    cancer_prop = cancer_prop
  ))
}


benchmark_reference <- function(deconv_ref, pseudobulk_list,
                                split_cancer = FALSE,
                                deconv_method = "qp") {
  deconv_res_list <- pseudobulk_list %>%
    lapply(
      deconvolute_pseudobulk,
      deconv_ref,
      split_cancer = split_cancer,
      method = deconv_method
    ) %>%
    extract(!unlist(lapply(., is.null)))

  deconv_prop_list <- deconv_res_list %>%
    lapply(function(deconv_res) {
      deconv_res$res
    })

  # Compute deconvolution errors
  deconv_err_vec <- deconv_prop_list %>%
    lapply(
      function(deconv_prop_df) {
        with(deconv_prop_df, rmse(prop_true, prop_deconv))
      }
    ) %>%
    unlist()

  print(mean(deconv_err_vec))

  # Generate overall result df
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

  # Summarize residuals
  deconv_residuals_df <- deconv_res_list %>%
    lapply(function(deconv_res) {
      deconv_res$residuals %>%
        {
          data.frame(
            transcript = rownames(.),
            residual = .
          )
        }
    }) %>%
    bind_rows(.id = "sample")

  resid_sum_df <- deconv_residuals_df %>%
    group_by(sample) %>%
    summarize(
      sum_sq_resid = sum(residual^2),
      sum_abs_resid = sum(abs(residual)),
      sum_resid = sum(residual)
    )

  # Compute cancer comp df
  cancer_prop_from_resid_df <- deconv_res_list %>%
    lapply(function(deconv_res) {
      deconv_res$cancer_prop
    }) %>%
    bind_rows(.id = "sample")

  cancer_prop_df <- all_prop_df %>%
    filter(celltype == "Cancer Epithelial") %>%
    select(sample, prop_true, prop_deconv, abs_err)

  cancer_comp_df <- cancer_prop_from_resid_df %>%
    left_join(cancer_prop_df, by = "sample") %>%
    left_join(resid_sum_df, by = "sample")

  return(list(
    "errors" = deconv_err_vec,
    "deconv_res" = all_prop_df,
    "deconv_sum" = all_prop_sum_df,
    "cancer_comp" = cancer_comp_df
  ))
}
