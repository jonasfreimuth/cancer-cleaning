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


norm_vec <- function(count_vec, type, base = 10, scale = 10^6) {
  vec_sum <- sum(count_vec)

  if (vec_sum == 0) {
    return(count_vec)
  }

  norm_vec <- count_vec / vec_sum

  # TODO Check if this check if log taking should be performed is not too
  # fragile.
  if (str_detect(type, "^log")) {
    norm_vec <- log(norm_vec + 1, base = base)
  }

  norm_vec <- norm_vec * scale

  return(norm_vec)
}

quantnorm_mat <- function(x, MAR = 2) {
  # FIXME This assumes a simple 2D matrix, if it gets more complicated this
  # needs to be adapted.
  CMAR <- 2

  if (MAR != 1) {
    CMAR <- 1
  }

  # Get per rank averages.
  rank_avg <- x %>%
    # Ensure rank order
    apply(CMAR, sort, na.last = TRUE) %>%
    t() %>%
    apply(MAR, mean, na.rm = TRUE) %>%
    # FIXME Check if this is necessary
    sort()

  # Match rank averages back to individual rows
  # In cases where multiple cells have the same ranks,
  # use the mean rank average for all tied cells.
  x_normed <- x %>%
    apply(MAR, function(vec, rank_avg) {
      avg_df <- data.frame(
        rank_avg = rank_avg,
        rank = rank(rank_avg)
      )

      # TEHE let's build a df for EVERY col in the data, I'm sure that's great
      # for performance.
      data.frame(
        rank = rank(vec, na.last = TRUE, ties.method = "first"),
        group = rank(vec, na.last = TRUE)
      ) %>%
        left_join(avg_df, by = "rank") %>%
        group_by(group) %>%
        summarize(
          group_avg = mean(rank_avg),
          group_ranks = paste(rank, collapse = "_"),
          .groups = "drop_last"
        ) %>%
        separate_rows(
          group_ranks,
          sep = "_"
        ) %>%
        extract2("group_avg") %>%
        return()
    },
    rank_avg = rank_avg
    )

  if (MAR == 1) {
    x_normed %<>%
      t()
  }

  dimnames(x_normed) <- dimnames(x)

  return(x_normed)
}


normalize_count_mat <- function(count_mat, type = "lognorm", ...) {
  # FIXME Be a bit more clever with matchig type to normalization performed.
  if (is.null(type)) {
    # NULL means do nothing.
    return(count_mat)
  }

  if (type %in% c("lognorm", "norm")) {
    norm_mat <- count_mat %>%
      apply(
        2,
        norm_vec,
        simplify = TRUE,
        type = type,
        ...
      )
  } else if (type == "tpm") {
    norm_mat <- count_mat %>%
      scuttle::calculateTPM(...)
  } else if (type == "logtpm") {
    warning(paste0(
      "The logtpm method is experimental. Double check whether it does what ",
      "you want."
    ))
    norm_mat <- count_mat %>%
      scuttle::calculateTPM(...) %>%
      # TODO Check if this should not be integrated into the tpm fun.
      add(1) %>%
      apply(2, log2)
  } else if (type == "quantile") {
    norm_mat <- count_mat %>%
      # as.matrix() %>%
      quantnorm_mat()
  } else {
    stop(paste0(
      "Normalization type \"", type, "\" is not supported."
    ))
  }

  return(norm_mat)
}
