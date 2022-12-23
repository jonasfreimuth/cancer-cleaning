library("here")
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
        scale_factor = 1
      )
    ),
    initialize = function(count_matrix, cell_indices,
                          cancer_celltypes = NULL) {
      stopifnot(
        all(c("CountMatrix", "R6") %in% class(count_matrix))
      )
      private$.count_matrix <- count_matrix
      private$.cell_indices <- cell_indices

      self$params$cancer_celltypes <- cancer_celltypes

      cancer_indices <- which(
        count_matrix$celltypes %in% cancer_celltypes
      )
      private$.cancer_indices <- cancer_indices
      private$.clean_indices <- setdiff(
        cell_indices,
        cancer_indices
      )
    }
  ),
  active = list(
    matrix_raw = function() {
      private$.count_matrix$matrix
    },
    matrix_clean = function() {
      # TODO Change this to lazy field.
      self$matrix_raw %>%
        magrittr::extract(
          i = ,
          j = private$.clean_indices,
          drop = FALSE
        )
    },
    matrix_cancer = function() {
      # TODO Change this to lazy field.
      self$matrix_raw  %>%
        magrittr::extract(
          i = ,
          j = private$.cancer_indices,
          drop = FALSE
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
    }
  ),
  private = list(
    .count_matrix = NULL,
    .cell_indices = NULL,
    .cancer_indices = NULL,
    .clean_indices = NULL
  )
)
