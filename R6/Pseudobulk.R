ibrary("here")
library("R6")

here::i_am("R6/Pseudobulk.R")

source(here("R6/CountMatrix.R"))

Pseudobulk <- R6Class(
  "Pseudobulk",
  public = list(
    params = list(
      id = NULL,
      cancer_celltypes = NULL,
      normalization = list(
        type = NULL,
        scale_factor = NULL
      )
    ),
    initialize = function(count_matrix, cell_indices,
                          cancer_celltypes = NULL) {
      stopifnot(
        all(c("CountMatrix", "R6") %in% class(count_matrix))
      )
      private$count_matrix <- count_matrix
      private$cell_indices <- cell_indices

      self$params$cancer_celltypes <- cancer_celltypes
    }
  ),
  active = list(
    matrix_raw = function() {
      private$count_marix$matrix
    },
    matrix = function() {
      self$matrix_raw %>%
        extract(
          i = ,
          j = private$cell_indices,
          drop = FALSE
        )
    },
    meta_raw = function() {
      private$count_matrix$meta
    },
    meta = function() {
      self$meta_raw %>%
        filter(cell %in% colnames(self$count_matrix))
    },
    transcripts = function() {
      rownames(self$matrix)
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
      self$matrix %>%
        rowSums() %>%
        norm_vec(
          type = self$normalization$type,
          scale = self$normalization$scale_factor
        )
    },
    transcript_abundances_cancer = function() {
      cancer_idcs <- self$meta$celltypes %in% self$params$cancer_celltypes

      self$matrix[, cancer_idcs] %>%
        rowSums() %>%
        norm_vec(
          type = self$normalization$type,
          scale = self$normalization$scale_factor
        )
    },
    celltype_counts = function() {
      self$meta$celltypes %>%
        table() %>%
        {
          # Conversion from table to named vector.
          set_names(as.vector(.), names(.))
        }
    }
  ),
  private = list(
    count_matrix = NULL,
    cell_indices = NULL
  )
)
