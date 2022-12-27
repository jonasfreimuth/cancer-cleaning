library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Pseudobulk.R")
)

source(here("R6/CountMatrix.R"))
source(here("R6/CountMatrixWrapper.R"))

source(here("functions/norm_functions.R"))

Pseudobulk <- R6Class(
  "Pseudobulk",
  inherit = CountMatrixWrapper,
  public = list(
    initialize = function(count_matrix, params) {
      stopifnot(
        all(c("CountMatrix", "R6") %in% class(count_matrix))
      )
      private$.check_params(params)
      private$.check_params_count_matrix_fit(params, count_matrix)

      super$initialize(count_matrix)

      private$.params <- params

      cancer_indices <- which(
        count_matrix$celltypes %in% params$cancer_celltypes
      )
      private$.cancer_indices <- cancer_indices
      private$.clean_indices <- setdiff(
        params$cell_indices,
        cancer_indices
      )
    }
  ),
  active = list(
    params = function() {
      private$.params
    },
    # TODO Fix issues with different matrix versions:
    # * Core Pseudobulk is per celltype transcript props normalized
    #    (with or without cancer cells?).
    # --> Celltype abundance matrix normalized (with or without cancer cells?).
    # *
    matrix_raw = function() {
      # TODO Add normalization downstream
      private$.count_matrix$matrix_orig
    },
    matrix_clean = function() {
      # TODO Change this to lazy field.
      self$matrix_raw %>%
        magrittr::extract(
          i = ,
          j = private$.clean_indices,
          drop = FALSE
        ) %>%
        normalize_count_mat(
          type = self$params$normalization$type,
          scale = self$params$normalization$scale_factor
        )
    },
    matrix_cancer = function() {
      # TODO Change this to lazy field.
      self$matrix_raw %>%
        magrittr::extract(
          i = ,
          j = private$.cancer_indices,
          drop = FALSE
        ) %>%
        normalize_count_mat(
          type = self$params$normalization$type,
          scale = self$params$normalization$scale_factor
        )
    },
    meta_raw = function() {
      private$.count_matrix$meta
    },
    meta = function() {
      self$meta_raw %>%
        magrittr::extract(
          i = private$.clean_indices,
          j = ,
          drop = FALSE
        )
    },
    transcripts = function() {
      # TODO Should zero sum transcripts be removed?
      rownames(self$matrix_clean)
    },
    bulk = function() {
      self$transcript_abundances %>%
        set_names(self$transcripts)
    },
    df = function() {
      self$bulk %>%
        {
          data.frame(
            IDs = names(.),
            sample = .
          )
        }
    },
    transcript_abundances = function() {
      self$matrix_clean %>%
        rowSums() %>%
        norm_vec(
          type = self$params$normalization$type,
          scale = self$params$normalization$scale_factor
        )
    },
    transcript_abundances_cancer = function() {
      # TODO Check if this makes sense wrt normalization.
      self$matrix_cancer %>%
        rowSums() %>%
        norm_vec(
          type = self$params$normalization$type,
          scale = self$params$normalization$scale_factor
        )
    },
    celltype_counts = function() {
      self$meta$celltype %>%
        table() %>%
        {
          # Conversion from table to named vector.
          set_names(as.vector(.), names(.))
        }
    },
    celltype_abundances = function() {
      self$celltype_counts %>%
        norm_vec(
          type = "norm",
          scale = 1
        )
    },
    celltype_matrix_clean = function() {
      # TODO Implement.
      stop("Not implemented")
    }
  ),
  private = list(
    .count_matrix = NULL,
    .cancer_indices = NULL,
    .clean_indices = NULL,
    .params = NULL,
    .check_params = function(params) {
      stopifnot(
        all(c("R6", "Params", "PseudobulkParams") %in% class(params))
      )
    },
    .check_params_count_matrix_fit = function(params, count_matrix) {
      stopifnot(
        length(params$cell_indices) <= count_matrix$n_cells,
        max(params$cell_indices) <= count_matrix$n_cells
      )
    }
  )
)
