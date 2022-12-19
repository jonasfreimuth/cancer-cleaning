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
