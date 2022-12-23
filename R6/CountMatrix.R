library("here")
library("data.table")
library("magrittr")
library("dplyr")
library("R6")

here::i_am("R6/CountMatrix.R")

source(here("functions/norm_functions.R"))

CountMatrix <- R6Class(
  "CountMatrix",
  public = list(
    params = list(
      normalization = list(
        type = NULL,
        scale_factor = 1
      )
    ),
    initialize = function(matrix, meta, downsample,
                          downsample_frac) {
      meta <- meta %>%
        # Ordering ensured here.
        arrange(cell)

      # Ordering ensured here.
      matrix <- matrix[, order(colnames(matrix))]

      private$check_matrix(matrix)
      private$check_meta(meta)
      private$check_matrix_meta_fit(matrix, meta)

      private$matrix <- matrix
      private$meta <- meta

      if (is.null(downsample_frac)) {
        private$downsample(downsample_frac, downsample_frac)
      }
    },
    get_matrix = function() {
      private$normalize_mat()
    },
    get_matrix_orig = function() {
      private$matrix
    }
  ),
  active = list(
    cells = function() {
      colnames(private$matrix)
    },
    celltypes = function() {
      private$meta$celltypes
    },
    transcripts = function() {
      rownames(private$matrix)
    },
    n_cells = function() {
      self$cells %>%
        # NOTE Cells should always be unique, so there is no safeguard here.
        length()
    },
    n_celltypes = function() {
      self$celltypes %>%
        unique() %>%
        length()
    },
    celltype_count_matrix = function() {
      private$matrix %>%
        t() %>%
        as.matrix() %>%
        rowsum(group = self$celltypes) %>%
        t() %>%
        normalize_count_mat(
          type = self$params$normalization$type,
          scale = self$params$normalization$scale_factor
        )
    }
  ),
  private = list(
    matrix = NULL,
    # TODO Ensure meta is sorted by matrix order of cells
    meta = NULL,
    normalize_mat = function() {
      normalize_count_mat(
        count_mat = private$matrix,
        type = self$params$normalization$type,
        scale = self$params$normalization$scale_factor
      )
    },
    downsample = function(cell_frac, transcript_frac) {
      rnd_cell_idx <- self$n_cells %>%
        {
          sample(x = seq_len(.), size = cell_frac * .)
        }
      rnd_transcript_idx <- self$n_celltypes %>%
        {
          sample(x = seq_len(.), size = transcript_frac * .)
        }

      private$matrix <- private$matrix[rnd_cell_idx, rnd_transcript_idx]

      private$meta <- private$meta %>%
        filter(cell %in% private$cells)
    },
    check_matrix = function(matrix) {
      stopifnot(
        nrow(matrix) > 0,
        ncol(matrix) > 0,
        length(colnames(matrix)) > 0,
        length(rownames(matrix)) > 0,
        is(matrix, "sparseMatrix")
      )
    },
    check_meta = function(meta) {
      stopifnot(
        is.data.frame(meta),
        c("cell", "celltype") %in% names(meta)
      )
    },
    check_matrix_meta_fit = function(matrix, meta) {
      stopifnot(
        ncol(matrix) == nrow(meta),
        all(colnames(matrix) == meta$cells)
      )
    }
  )
)
