# TODO Clean up loaded pacakges.
library("purrr")
library("stringr")
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


clean_sigmat <- function(sigmat) {
  # TODO Include uniquely uncounted transcripts
  .row_is_identifiying <- function(row, frac = 0) {
    # Determine if a row can be used to identify a column. Whether or not that
    # is the case is determined by whether a fraction of elements is either
    # non-zero or zero. The fraction maps to the number of elements which is at
    # minimum 1 and at maximum the length of the row.
    .test_row <- function(row, frac) {
      len <- length(row)
      len_allowed <- len * frac

      if (len_allowed < 1) {
        len_allowed <- 1
      } else if (len_allowed > len) {
        len_allowed <- len
      }

      row_sum <- sum(row > 0)

      row_sum <= len_allowed && row_sum > 0
    }

    test_row_list <- list(
      f = row,
      r = 1 - row
    )

    res <- test_row_list %>%
      lapply(
        .test_row,
        frac = frac
      ) %>%
      purrr::reduce(or)

    return(res)
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


reference_from_thresh <- function(count_thresh,
                                  proto_sigmat) {
  sigmat <- proto_sigmat %>%
    is_greater_than(count_thresh) %>%
    # Simple as.numeric() returns a vector.
    multiply_by(1) %>%
    clean_sigmat() %>%
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


hampel_intervall <- function(x, mad_mult = 3) {
  x_median <- median(x)
  x_mad_mult <- mad(x) * mad_mult

  return(c(x_median - x_mad_mult, x_median + x_mad_mult))
}


is_uniform <- function(x) {
  all(x == x[1])
}


clean_nbin_sigmat <- function(sigmat, trim_used = 0.1) {
  ctypes <- colnames(sigmat)
  # remove rows that are all the same
  sigmat_wip <- sigmat %>%
    extract(i = !apply(., 1, is_uniform), j = , drop = FALSE)

  # calc per col/ctype outliers
  outlier_intervals <- sigmat_wip %>%
    apply(
      2,
      hampel_intervall
    )

  outlying_transcripts <- ctypes %>%
    # assign each col its range
    lapply(
      function(celltype, sigmat, interval_mat) {
        sel_vec <- colnames(sigmat) == celltype
        list(
          count_vec = sigmat[, sel_vec],
          range_vec = interval_mat[, sel_vec]
        )
      },
      sigmat_wip,
      outlier_intervals
    ) %>%
    set_names(ctypes) %>%
    # extract the trim fraction of most extreme outliers
    lapply(
      function(count_range_el, trim) {
        count_vec <- count_range_el %>%
          extract2("count_vec")
        range_vec <- count_range_el %>%
          extract2("range_vec")

        outlier_list <- count_vec %>%
          {
            list(
              "lower" = extract(., . < range_vec[1]) %>%
                sort(),
              "upper" = extract(., . > range_vec[2]) %>%
                sort(decreasing = TRUE)
            )
          } %>%
          lapply(
            function(outlier_vec, trim) {
              trim_cut <- length(outlier_vec) * (trim / 2)
              outlier_vec[seq_len(trim_cut)]
            },
            trim
          )

        outlier_vec <- unlist(outlier_list) %>%
          set_names(str_replace_all(
            names(.),
            "(upper.|lower.)",
            ""
          ))

        return(outlier_vec)
      },
      trim_used
    ) %>%
    lapply(names) %>%
    unlist() %>%
    unique()

  sigmat_wip <- sigmat_wip %>%
    extract(i = rownames(.) %in% outlying_transcripts, j = , drop = FALSE)

  # use top / botton n percent?
  sigmat_clean <- sigmat_wip
  return(sigmat_clean)
}


reference_non_binary <- function(proto_sigmat) {
  sigmat <- proto_sigmat %>%
    clean_nbin_sigmat() %>%
    dedupe_sigmut_mat()

  # sigmats with colSums equal 0 lead to deconv troubles
  # TODO check what problems this causes, e.g. leading to
  # 1 or 0 col sigmats.
  col_sums <- colSums(sigmat)

  sigmat <- sigmat[, col_sums > 0, drop = FALSE]

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

  # NOTE This relies on celltype_map being ordered by count_mat cells
  cancer_idx <- intersect(
    idx_vec,
    which(names(celltype_map) == "Cancer Epithelial")
  )

  cancer_expression <- count_mat[, cancer_idx] %>%
    rowSums()

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
      celltype_counts = celltype_counts,
      transcript_counts_cancer = cancer_expression
    )
  )
}


split_deconv_res <- function(deconv_prop_vec, sep = "_") {
  if (is.null(names(deconv_prop_vec))) {
    stop(paste0(
      "deconv_prop_vec needs to be a named numeric vector, but is unnamed."
    ))
  }

  if (!any(str_detect(names(deconv_prop_vec), sep))) {
    return(deconv_prop_vec)
  }

  prop_df <- data.frame(
    name = names(deconv_prop_vec),
    prop = deconv_prop_vec
  ) %>%
    mutate(group_size = str_count(name, sep) + 1) %>%
    separate_rows(name, sep = sep) %>%
    mutate(prop = prop / group_size)

  deconv_prop_split <- prop_df %>%
    extract2("prop") %>%
    set_names(prop_df$name)

  return(deconv_prop_split)
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
  # Subsetting is necessary as the residuals' length won't match otherwise.
  transcript_props_cancer <- pseudobulk[["transcript_counts_cancer"]] %>%
    extract(names(.) %in% deconv_ref$IDs) %>%
    divide_by(sum(.)) %>%
    extract(!is.na(.))

  celltype_counts <- pseudobulk[["celltype_counts"]]

  if (split_cancer) {
    celltype_counts <- celltype_counts %>%
      extract(!(names(.) %in% cancer_cols))
  }

  celltype_props <- celltype_counts %>%
    divide_by(sum(.)) %>%
    extract(!is.na(.))

  input_list <- list(
    tp = transcript_props,
    cp = celltype_props,
    kp = transcript_props_cancer
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

  if (split_cancer) {
    ref_names_clean <- deconv_ref %>%
      names() %>%
      str_replace_all(
        pattern = paste0(
          "(", paste(cancer_cols, collapse = "|"), ")"
        ),
        replacement = ""
      ) %>%
      str_replace_all(
        pattern = "(^_|(?<=_)_|_$)",
        replacement = ""
      )

    clean_ref_name_idx <- ref_names_clean != ""

    deconv_ref <- deconv_ref %>%
      extract(clean_ref_name_idx) %>%
      set_names(ref_names_clean[clean_ref_name_idx])
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

  cancer_expr_cor <- cor(transcript_props_cancer, as.vector(deconv_resid))

  resid_expr_df <- data.frame(
    transcript = names(transcript_props_cancer),
    cancer_expr = transcript_props_cancer,
    resid = deconv_resid,
    deconv_pred = transcript_props_pred
  )

  true_prop_df <- celltype_props %>%
    {
      data.frame(celltype = names(.), prop = .)
    } %>%
    separate_rows(celltype, sep = "_")

  deconv_prop_df <- dp_vec %>%
    split_deconv_res() %>%
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
    cancer_expr_corr = cancer_expr_cor,
    resid_expr_df = resid_expr_df
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
    lapply(
      extract2,
      "res"
    )

  # Compute deconvolution errors
  deconv_err_vec <- deconv_prop_list %>%
    lapply(
      function(deconv_prop_df) {
        with(deconv_prop_df, rmse(prop_true, prop_deconv))
      }
    ) %>%
    unlist()

  deconv_corr_df <- deconv_res_list %>%
    lapply(
      extract2,
      "cancer_expr_corr"
    ) %>%
    unlist() %>%
    {
      data.frame(
        sample = as.character(seq_along(.)),
        cancer_expr_corr = .
      )
    }

  resid_expr_df <- deconv_res_list %>%
    dfextract("resid_expr_df", "sample")

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
      unq_n_ct = paste(unique(n_ctypes), sep = "_")
    ) %>%
    left_join(deconv_corr_df, by = "sample") %>%
    mutate(rmse = deconv_err_vec)

  return(list(
    "deconv_res" = all_prop_df,
    "deconv_sum" = all_prop_sum_df,
    "resid_expr_df" = resid_expr_df
  ))
}
