ibrary("here")
library("R6")

here::i_am("R6/Pseudobulk.R")

source(here("R6/CountMatrix"))


Pseudobulk <- R6Class(
  "Pseudobulk",
  inherit = CountMatrix,
  public = list(
    params = list(
      id = NULL,
      # TODO Change this from logical to char vec giving cols to split.
      split_cancer = NULL,
      cancer_celltypes = NULL,
      normalization = list(
        type = NULL,
        scale_factor = NULL
      )
    ),
    initialize = function(experiment, cell_indices,
                          split_cancer,
                          cancer_celltypes,
                          normalization_type = "lognorm",
                          scale_factor = 1) {
      # TODO Check experiment type.
      # TODO Check cexperimetn matrix and meta ordering.
      self$matrix <- experiment$matrix[, cell_indices]
      self$meta <- experiment$meta %>%
        filter(cell %in% self$cells)

      self$params$split_cancer <- split_cancer
      self$params$cancer_celltypes <- cancer_celltypes

      self$params$normalization$type <- normalization_type
      self$params$normalization$scale_factor <- scale_factor
    }
  ),
  active = list(
    transcript_counts = function() {
      self$matrix %>%
        rowSums() %>%
        norm_vec(
          type = norm_type,
          scale = scale
        )
    },
    transcript_counts_cancer = function() {
      cancer_idcs <- count_matrix$celltypes %in% self$params$cancer_celltypes

      self$matrix[, cancer_idcs] %>%
        rowSums() %>%
        norm_vec(
          type = norm_type,
          scale = scale
        )
    },
    celltype_counts = function() {
      self$celltypes %>%
        table() %>%
        {
          # Conversion from table to named vector.
          set_names(as.vector(.), names(.))
        }
    }
  )
)
