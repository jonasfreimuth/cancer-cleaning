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

# TODO Once it is clear which functions are actually useful, add proper
# docstrings.
dfextract <- function(list, name, .id) {
  list %>%
    lapply(
      function(list_el) {
        list_el %>%
          extract2(name)
      }
    ) %>%
    bind_rows(.id = .id)
}


dfextract2 <- function(list, name, .id_inner, .id_outer) {
  list %>%
    lapply(
      dfextract,
      name,
      .id = .id_inner
    ) %>%
    bind_rows(.id = .id_outer)
}


load_experiment <- function(count_mat_file, rowname_file, colname_file,
                            meta_file, testing = FALSE) {
  # FIXME Make this function more general
  data_full_meta <- fread(meta_file)

  if ("V1" %in% names(data_full_meta)) {
    data_full_meta %<>%
      dplyr::rename(cell = V1)
  }

  data_full_rownames <- readLines(rowname_file)

  data_full_colnames <- readLines(colname_file)

  data_full_matrix <- readMM(count_mat_file)

  data_full_matrix@Dimnames <- list(
    data_full_rownames,
    data_full_colnames
  )

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


clean_names <- function(names, cancer_names, sep = "_") {
  # Remove cancer_names substring from names, assuming concatenated names are
  # separe
  # Input:
  # * names: A vector of column names, potentially including cancer_names as
  #   substrings.
  # * cancer_names: A vector of strings (denoting cancer columns) that should
  #   be removed from names.
  # * sep: A string separating concatenated names. CAUTION: No sanitation is
  #   performed! (I can't be asked)
  # Output:
  # names, but with all occurrences of cancer_names removed. Also removed are
  # leading / trailing / duplicated '_' characters (The results of a removed
  # cancer_names substring), and empty entries removed.
  names %>%
    str_replace_all(
      pattern = paste0(
        "(", paste(cancer_names, collapse = "|"), ")"
      ),
      replacement = ""
    ) %>%
    # Handle collapsed sigmat cols from deduping.
    str_replace_all(
      pattern = str_glue("(^{sep}|(?<={sep}){sep}|{sep}$)"),
      replacement = ""
    ) %>%
    extract(. != "")
}


split_concat_celltype_props <- function(vec, sep = "_") {
  # Adapted from pigx_sars_cov_2 pipeline approach to splitting variant props
  # originating from variants with the same signature.
  .split_cprop_row <- function(prop_group, sep = "_") {
    group_members <- strsplit(prop_group["group_string"], sep)

    group_len <- length(group_members)

    # Assign each individual member the same fraction of the prop, for details
    # see the pigx_sars_cov_2 publication. Something about an assumption of
    # normal distribution among indistinguishable variants.
    member_props <- rep(
      as.numeric(prop_group["prop"]),
      group_len
    ) /
      group_len

    data.frame(
      group_string = group_members,
      prop = member_props
    )
  }

  prop_df <- vec %>%
    {
      data.frame(
        group_string = names(.),
        prop = .
      )
    } %>%
    # NOTE Inside the apply call we are dealing with a character matrix. Input
    # is still a dataframe because it is easier and we get colnames in the apply
    # call.
    # FIXME Make this more explicit by having the input be a character matrix
    # already.
    apply(
      1,
      .split_cprop_row,
      sep
    ) %>%
    bind_rows()

  prop_df$prop %>%
    set_names(prop_df$group_string)
}


frequency_bins <- function(x, prob_seq = seq(0, 1, 0.2)) {
  x %>%
    quantile(probs = prob_seq) %>%
    set_names(as.character(round(., 2)))
}

quantile_counts <- function(x, prob) {
  if (length(prob) > 2) {
    stop("prob needs to be either length 1 or 2.")
  }
  if (length(prob) == 1) {
    prob <- c(prob, 1 - prob)
  }

  bins <- quantile(x, prob)

  low_vec <- x <= bins[1]
  high_vec <- x >= bins[2]

  n_low <- sum(low_vec)
  n_mid <- sum(!low_vec & !high_vec)
  n_high <- sum(high_vec)

  c("low" = n_low, "mid" = n_mid, "high" = n_high)
}


is_quantile_unique <- function(x, prob = 0.2) {
  q_counts <- quantile_counts(x, prob)

  if (q_counts["mid"] > 0) {
    return(FALSE)
  }

  return(c(xor(q_counts["low"] == 1, q_counts["high"] == 1), use.names = FALSE))
}


mean_imbalance <- function(x) {
  mean_x <- mean(x)

  low <- x < mean_x
  high <- x > mean_x

  if (sum(low) == 1) {
    return(x[low] - mean_x)
  } else if (sum(high) == 1) {
    return(x[high] - mean_x)
  } else {
    return(0)
  }
}


mean_dists <- function(x) {
  mean_x <- mean(x)

  x - mean_x
}


hampel_interval <- function(x, mad_mult = 3) {
  x_median <- median(x)
  x_mad_mult <- mad(x) * mad_mult

  return(c(x_median - x_mad_mult, x_median + x_mad_mult))
}


outlier_dist <- function(x, mad_mult = 3) {
  dists <- mean_dists(x)
  dist_interval <- hampel_interval(dists, mad_mult)

  outlier_idcs <- dists < dist_interval[1] | dists > dist_interval[2]

  if (sum(outlier_idcs) == 1) {
    return(dists[outlier_idcs])
  } else {
    return(0)
  }
}


otsu_crit <- function(x, thresh) {
  x_thresh <- x >= thresh

  true_frac <- sum(x_thresh) / length(x_thresh)

  if (true_frac == 0 || true_frac == 1) {
    # Var won't be computable for low or high -> This threshold can't be
    # considered
    return(Inf)
  }

  x_high_var <- var(x[x_thresh])
  x_low_var <- var(x[!x_thresh])

  x_high_var * true_frac + x_low_var * (1 - true_frac)
}


otsu_thresh <- function(x) {
  # This function probably won't help in finding uniquely identifying
  # transcripts, due to the most extreme cases not working as
  # true_frac %in% c(0,1) in those cases.
  thresh_vec <- unique(x)

  crit_vec <- lapply(
    thresh_vec,
    otsu_crit,
    x = x
  ) %>%
    unlist()

  thresh_vec[ctrit_vec == min(crit_vec)]
}


seq_power <- function(start, stop, step_frac, power = 10) {
  seq_span <- stop - start

  if (seq_span == 0) {
    proto_seq <- 0
  } else {
    proto_seq <- seq(0, 1, step_frac) %>%
      raise_to_power(power) %>%
      {
        (. - min(.)) / (max(.) - min(.))
      }
  }
  return((proto_seq * seq_span) + start)
}
