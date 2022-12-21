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
    initialize = function(matrix, meta) {
      private$matrix <- matrix
      private$meta <- meta
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
      cell_colname <- "cell"
      celltype_colname <- "celltype_major"

      celltype_map <- private$meta %>%
        extract(c(cell_colname, celltype_colname))

      dataframe(cells = self$cells) %>%
        left_join(celltype_map, by = c("cells", cell_colname)) %>%
        extract2(celltype_colname)
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
    }
  )
)
