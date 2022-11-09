library("dplyr")
library("magrittr")

dedupe_sigmut_mat <- function(sigmut_mat, var_sep = "_") {
  # input:
  #   * sigmut_mat:
  #     A matrix with rows corresponing to singnature mutations and cols
  #     corresponding to variants, each dimension named. Entries are 1 if the
  #     mutation of the row is a signature mutation of the variant in the col,
  #     and 0 otherwise.
  #   * var_sep:
  #     A string separating the names of variables with identical columns.
  # output:
  #   * The input matrix with identical columns and their colnames merged.
  #   * A list of character vectors, giving names of variants which have equal
  #     cols.

  # Handle case of only one col present (might get read as a vector)
  if (is.null(dim(sigmut_mat))) {
    return(sigmut_mat)
  }

  variant_names <- colnames(sigmut_mat)

  is_dupe <- duplicated(sigmut_mat, MARGIN = 2)

  if (any(is_dupe)) {
    dupe_variants <- variant_names[is_dupe]

    # coerce back to dataframe for easier processing
    sigmut_mat_df <- as.data.frame(sigmut_mat)

    # find out of which variant a dupe variant is a dupe of, generate groups
    # of variants which are duplicates of each other
    dupe_group_list <- list()

    for (dupe_var in dupe_variants) {
      if (!dupe_var %in% unique(unlist(dupe_group_list))) {
        dupe_var_col <- sigmut_mat_df[[dupe_var]]

        dupe_group_logi <- apply(
          sigmut_mat,
          2,
          function(col, dupe_col) {
            all(col == dupe_col)
          }, dupe_var_col
        )

        dupe_group_vec <- variant_names[dupe_group_logi]

        dupe_group_list[[dupe_group_vec[1]]] <- dupe_group_vec
      }
    }

    # concat dupe groups to form a new composite name for the now unique col
    dupe_group_names <- lapply(dupe_group_list, paste, collapse = var_sep) %>%
      unlist()

    # juggle names to get a named vector with names and values flipped
    # needed by dplyr::rename()
    old_names <- names(dupe_group_names)
    new_names <- dupe_group_names

    dupe_group_names <- old_names %>%
      set_names(new_names)

    # generate deduped signature matrix
    # is a col was duplicated this contains only the first col of each dupe
    # group
    dedupe_res <- sigmut_mat_df[, !is_dupe, drop = FALSE] %>%
      rename(!!dupe_group_names) %>%
      replace(is.na(.), 0)
  } else {
    dedupe_res <- sigmut_mat %>%
      as.data.frame()
  }

  return(dedupe_res)
}
