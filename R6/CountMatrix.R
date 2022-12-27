library("here")
library("data.table")
library("magrittr")
library("dplyr")
library("R6")

suppressMessages(
  here::i_am("R6/CountMatrix.R")
)

source(here("functions/norm_functions.R"))

CountMatrix <- R6Class(
  "CountMatrix",
  public = list(
    initialize = function(matrix, meta, downsample_frac = NULL) {
      meta <- meta %>%
        # Ordering ensured here.
        arrange(cell)

      # Ordering ensured here.
      matrix <- matrix[order(rownames(matrix)), order(colnames(matrix))]

      private$.check_matrix(matrix)
      private$.check_meta(meta)
      private$.check_matrix_meta_fit(matrix, meta)

      private$.matrix <- matrix
      private$.meta <- meta

      if (!is.null(downsample_frac)) {
        private$.downsample(downsample_frac, downsample_frac)
      }
    }
  ),
  active = list(
    params = function() {
      private$.params
    },
    matrix_orig = function() {
      private$.matrix
    },
    meta = function() {
      private$.meta
    },
    cells = function() {
      colnames(private$.matrix)
    },
    celltypes = function() {
      private$.meta$celltype
    },
    transcripts = function() {
      rownames(private$.matrix)
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
    n_transcripts = function() {
      self$transcripts %>%
        length()
    },
    celltype_count_matrix = function() {
      private$.matrix %>%
        t() %>%
        as.matrix() %>%
        rowsum(group = self$celltypes) %>%
        t()
    }
  ),
  private = list(
    .matrix = NULL,
    # TODO Ensure meta is sorted by matrix order of cells
    .meta = NULL,
    .downsample = function(cell_frac, transcript_frac) {
      rnd_cell_idx <- self$n_cells %>%
        {
          sample(x = seq_len(.), size = cell_frac * .)
        }
      rnd_transcript_idx <- self$n_transcripts %>%
        {
          sample(x = seq_len(.), size = transcript_frac * .)
        }

      private$.matrix <- private$.matrix[rnd_transcript_idx, rnd_cell_idx]

      private$.meta <- private$.meta %>%
        filter(cell %in% colnames(private$.matrix))
    },
    .check_matrix = function(matrix) {
      stopifnot(
        nrow(matrix) > 0,
        ncol(matrix) > 0,
        length(colnames(matrix)) > 0,
        length(rownames(matrix)) > 0,
        is(matrix, "sparseMatrix")
      )
    },
    .check_meta = function(meta) {
      stopifnot(
        is.data.frame(meta),
        c("cell", "celltype") %in% names(meta)
      )
    },
    .check_matrix_meta_fit = function(matrix, meta) {
      stopifnot(
        ncol(matrix) == nrow(meta),
        all(colnames(matrix) == meta$cells)
      )
    }
  )
)
