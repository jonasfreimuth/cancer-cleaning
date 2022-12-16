# TODO Clean up loaded pacakges.
library("here")
library("glmGamPoi")
library("DESeq2")
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

source(here("functions/norm_functions.R"))


get_de_transcripts <- function(count_mat, meta, design) {
  # TODO Get the functions in here to STFU.
  meta <- meta %>%
    arrange(cell) %>%
    # Prevent DESeq fun from complaining
    mutate(across(where(is.character), as.factor))
  count_mat <- count_mat[, order(colnames(count_mat))]

  # TODO Is prescaling useful here?
  size_factors <- count_mat %>%
    scuttle::pooledSizeFactors()

  ds2_data <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = meta,
    design = design
  )

  if (length(design) >= 3) {
    # This assumes the last term of the design formula is the only one of
    # interest.
    # TODO Enforce this assumption.
    reduced_design <- design[-length(design)]
  } else {
    reduced_design <- ~1
  }

  ds2_data %>%
    `sizeFactors<-`(value = size_factors) %>%
    DESeq(
      # Args are set according to
      # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis.
      test = "LRT",
      fitType = "glmGamPoi",

      # This is recommended, but according to the source above shouldn't do
      # anything with test = "LRT".
      useT = TRUE,

      # This should also be already set when using fitType = "glmGamPoi".
      minmu = 10^-6,
      reduced = reduced_design,
      quiet = TRUE
    ) %>%
    results(
      tidy = TRUE
    ) %>%
    # Caution, results may contain NAs, seems to mainly depend on whether
    # zero cols were cleaned beforehand.
    return()
}


count_df_from_reference <- function(reference, feature_scale = FALSE) {
  if (!"IDs" %in% names(reference)) {
    stop(paste0(
      "IDs col needs to be present in reference, as it is assumed to contain ",
      "transcript names."
    ))
  }

  # Input must be a reference df, IDs needs to be present
  count_df <- reference %>%
    set_rownames(.$IDs) %>%
    select(-IDs) %>%
    t() %>%
    as.data.frame() %>%
    mutate(
      celltype = rownames(.)
    ) %>%
    select(celltype, everything()) %>%
    pivot_longer(
      !all_of(c("celltype")),
      names_to = "transcript",
      values_to = "abundance"
    )

  if (feature_scale) {
    count_df <- count_df %>%
      group_by(transcript) %>%
      mutate(abundance = feature_scale(abundance))
  }

  return(count_df)
}


sigmat_qc_plot <- function(reference, title = NULL, feature_scale = FALSE) {
  count_df <- count_df_from_reference(reference, feature_scale)

  celltype_means <- count_df %>%
    group_by(celltype) %>%
    summarize(
      mean = mean(abundance),
      .groups = "drop_last"
    ) %>%
    # Align to count_df for reorder
    {
      left_join(select(count_df, celltype), ., by = "celltype")
    } %>%
    extract2("mean")

  transcript_means <- count_df %>%
    group_by(transcript) %>%
    summarize(
      mean = mean(abundance),
      .groups = "drop_last"
    ) %>%
    # Align to count_df for reorder
    {
      left_join(select(count_df, transcript), ., by = "transcript")
    } %>%
    extract2("mean")


  qc_plot <- ggplot(
    count_df,
    aes(
      x = reorder(celltype, -celltype_means),
      y = reorder(transcript, transcript_means), fill = abundance
    )
  ) +
    geom_raster() +
    labs(
      x = "Celltype",
      y = "Transcript",
      title = title
    ) +
    scale_fill_gradient(
      name = "Abundance",
      low = "dodgerblue4",
      high = "yellow"
    ) +
    theme_benchmark() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_blank(),
      legend.text = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom"
    )

  return(qc_plot)
}


outlier_plot <- function(reference, selection) {
  reference %<>%
    filter(IDs %in% selection)

  count_df <- count_df_from_reference(reference, feature_scale = TRUE)

  outlier_plot <- ggplot(count_df, aes(transcript, abundance)) +
    geom_boxplot() +
    labs(
      x = "Transcript",
      y = "Per celltype abundance"
    ) +
    theme_benchmark()
}


remove_unidentifying_bin_rows <- function(sigmat) {
  .row_is_identifying <- function(row, frac = 0) {
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
      .row_is_identifying
    )

  sigmat_unique <- sigmat[unique_transcript_idx_vec, , drop = FALSE]

  return(sigmat_unique)
}


bin_reference_from_thresh <- function(count_thresh,
                                      proto_sigmat) {
  # TODO: Consider transcript counts as weights.
  # TODO: Add Others col?
  sigmat <- proto_sigmat %>%
    is_greater_than(count_thresh) %>%
    # Simple as.numeric() returns a vector.
    multiply_by(1) %>%
    remove_unidentifying_bin_rows() %>%
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


is_uniform <- function(x) {
  all(x == x[1])
}


reference_from_marker_df <- function(proto_sigmat, marker_df, thresh) {
  if (thresh == 0) {
    # Set it to a small number, otherwise the whole proto sigmat will be
    # returned.
    thresh <- 10^-16
  }
  marker_df_thresh <- marker_df %>%
    slice_max(metric, prop = thresh)

  # TODO Think about whether further cleaning is needed, maybe look at
  # removed Hampel stuff for that.
  proto_sigmat %>%
    extract(
      i = rownames(.) %in% marker_df_thresh$transcript, j = , drop = FALSE
    ) %>%
    as.data.frame() %>%
    mutate(IDs = rownames(.)) %>%
    select(IDs, everything()) %>%
    return()
}


create_celltype_map <- function(celltype, meta_df, cell_colname,
                                celltype_colname) {
  meta_df %>%
    filter(if_any(all_of(celltype_colname), ~ .x == celltype)) %>%
    as.data.frame() %>%
    extract2(cell_colname)
}


pseudobulk_from_idx <- function(idx_vec, count_mat, celltype_map,
                                norm_type = "tpm",
                                scale = norm_scale,
                                cancer_cols = "Cancer Epithelial") {
  bulk_count_mat <- count_mat[, idx_vec]

  transcript_counts <- bulk_count_mat %>%
    rowSums() %>%
    norm_vec(
      type = norm_type,
      scale = scale
    )

  # NOTE This relies on celltype_map being ordered by count_mat cells
  cancer_idx <- intersect(
    idx_vec,
    which(names(celltype_map) %in% cancer_cols)
  )

  cancer_expression <- count_mat[, cancer_idx] %>%
    rowSums() %>%
    norm_vec(
      type = norm_type,
      scale = scale
    )

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

  bulk_proto_sigmat <- bulk_count_mat %>%
    t() %>%
    as.matrix() %>%
    # NOTE This relies on celltype_map being ordered by count_mat cells
    rowsum(group = names(celltype_map)[idx_vec]) %>%
    t() %>%
    # TODO Check if this is right.
    normalize_count_mat(
      type = norm_type,
      scale = scale
    )

  return(
    list(
      transcript_counts = transcript_counts,
      celltype_counts = celltype_counts,
      transcript_counts_cancer = cancer_expression,
      bulk_proto_sigmat = bulk_proto_sigmat
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

  transcript_props_all <- pseudobulk[["transcript_counts"]]
  transcript_props_all_cancer <- pseudobulk[["transcript_counts_cancer"]]

  transcript_props_marker <- transcript_props_all %>%
    extract(names(.) %in% deconv_ref$IDs)
  transcript_props_marker_cancer <- transcript_props_all_cancer %>%
    extract(names(.) %in% deconv_ref$IDs)

  bulk_proto_sigmat <- pseudobulk[["bulk_proto_sigmat"]]

  celltype_counts <- pseudobulk[["celltype_counts"]]

  if (split_cancer) {
    celltype_counts <- celltype_counts %>%
      extract(!(names(.) %in% cancer_cols))
  }

  celltype_props <- celltype_counts %>%
    # NOTE The celltype_props need to be on the same scale as what we get from
    # deconvolution, that is why this normalization happens here.
    divide_by(sum(.)) %>%
    extract(!is.na(.))

  input_list <- list(
    cp = celltype_props
  )

  input_lens <- lapply(input_list, length)

  if (any(input_lens == 0)) {
    return(NULL)
  }

  deconv_bulk <- transcript_props_marker %>%
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

    # NOTE The entries corresponding to IDs col need to be dropped. They are
    # assumed to be at pos 1, hence the [-1] subsetting.
    bulk_proto_sigmat <- bulk_proto_sigmat %>%
      extract(i = , j = clean_ref_name_idx[-1], drop = FALSE) %>%
      set_colnames(ref_names_clean[-1][clean_ref_name_idx[-1]])
  }

  if (ncol(deconv_ref) <= 2) {
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

  # Marker residuals
  transcript_props_marker_pred <- sigmat %*% dp_vec

  deconv_resid_marker <- transcript_props_marker - transcript_props_marker_pred

  cor_marker_df <- data.frame(
    cancer_expr_v_resid = cor(
      transcript_props_marker_cancer, as.vector(deconv_resid_marker)
    ),
    cancer_expr_v_deconv_pred = cor(
      transcript_props_marker_cancer, as.vector(transcript_props_marker_pred)
    )
  )

  resid_expr_marker_df <- data.frame(
    transcript = names(transcript_props_marker_cancer),
    cancer_expr = transcript_props_marker_cancer,
    resid = deconv_resid_marker,
    deconv_pred = transcript_props_marker_pred
  )

  # All transcript residuals
  transcript_props_all_pred <- bulk_proto_sigmat %*% dp_vec

  deconv_resid_all <- transcript_props_all - transcript_props_all_pred

  cor_all_df <- data.frame(
    cancer_expr_v_resid = cor(
      transcript_props_all_cancer, as.vector(deconv_resid_all)
    ),
    cancer_expr_v_deconv_pred = cor(
      transcript_props_all_cancer, as.vector(transcript_props_all_pred)
    )
  )

  resid_expr_all_df <- data.frame(
    transcript = names(transcript_props_all_cancer),
    cancer_expr = transcript_props_all_cancer,
    resid = deconv_resid_all,
    deconv_pred = transcript_props_all_pred
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
    cor_marker_df = cor_marker_df,
    resid_expr_marker_df = resid_expr_marker_df,
    cor_all_df = cor_all_df,
    resid_expr_all_df = resid_expr_all_df
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


  deconv_corr_marker_df <- deconv_res_list %>%
    lapply(
      extract2,
      "cor_marker_df"
    ) %>%
    bind_rows() %>%
    mutate(
      sample = as.character(seq_len(nrow(.)))
    )

  deconv_corr_all_df <- deconv_res_list %>%
    lapply(
      extract2,
      "cor_all_df"
    ) %>%
    bind_rows() %>%
    mutate(
      sample = as.character(seq_len(nrow(.)))
    )


  resid_expr_marker_df <- deconv_res_list %>%
    dfextract("resid_expr_marker_df", "sample")

  resid_expr_all_df <- deconv_res_list %>%
    dfextract("resid_expr_all_df", "sample")


  cor_both_df <- left_join(
    deconv_corr_marker_df,
    deconv_corr_all_df,
    by = "sample",
    suffix = c("_marker", "_all")
  )


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
    left_join(cor_both_df, by = "sample") %>%
    mutate(rmse = deconv_err_vec)

  return(list(
    "deconv_res" = all_prop_df,
    "deconv_sum" = all_prop_sum_df,
    "resid_expr_marker_df" = resid_expr_marker_df,
    "resid_expr_all_df" = resid_expr_all_df
  ))
}
