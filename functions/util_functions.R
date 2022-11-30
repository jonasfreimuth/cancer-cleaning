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


seq_base <- function(start, stop, step_frac, base = 10) {
  if (is.null(base)) {
    return(seq(from = start, to = stop, by = step_frac * (stop - start)))
  }

  start <- start + 1
  stop <- stop + 1

  start_base <- log(start, base = base)
  stop_base <- log(stop, base = base)

  seq_len <- stop_base - start_base
  step_size <- seq_len * step_frac
  step_seq <- seq(from = start_base, to = stop_base, by = step_size)

  return(base^step_seq)
}
