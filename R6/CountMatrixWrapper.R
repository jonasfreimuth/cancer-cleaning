library("here")
library("R6")

suppressMessages(
  here::i_am("R6/CountMatrixWrapper.R")
)

source(here("R6/CountMatrix.R"))

# Class for facilitating access to core CountMatrix properties for inheriting
# objects.
CountMatrixWrapper <- R6Class(
  "CountMatrixWrapper",
  public = list(
    initialize = function(count_matrix) {
      private$.count_matrix <- count_matrix
    }
  ),
  active = list(
    matrix = function() {
      private$.count_matrix$matrix
    },
    matrix_orig = function() {
      private$.count_matrix$matrix_orig
    },
    meta = function() {
      private$.count_matrix$meta
    },
    cells = function() {
      private$.count_matrix$cells
    },
    celltypes = function() {
      private$.count_matrix$celltypes
    },
    transcripts = function() {
      private$.count_matrix$transcripts
    },
    n_cells = function() {
      private$.count_matrix$n_cells
    },
    n_celltypes = function() {
      private$.count_matrix$n_celltypes
    },
    n_transcripts = function() {
      private$.count_matrix$n_transcripts
    }
  ),
  private = list(
    .count_matrix = NULL,
    .check_count_matrix = function(count_matrix) {
      stopifnot(
        c("CountMatrix", "R6") %in% class(count_matrix)
      )
    }
  )
)
