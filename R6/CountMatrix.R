library("data.table")
library("magrittr")
library("dplyr")
library("R6")

CountMatrix <- R6Class(
  "CountMatrix",
  public = list(
    matrix = NULL,
    # TODO Ensure meta is sorted by matrix order of cells
    meta = NULL,
    initialize = function(matrix, meta) {
      self$matrix <- matrix
      self$meta <- meta
    },
  ),
  active = list(
    cells = function() {
      colnames(self$matrix)
    },
    celltypes = function() {
      cell_colname <- "cell"
      celltype_colname <- "celltype_major"

      celltype_map <- self$meta %>%
        extract(c(cell_colname, celltype_colname))

      dataframe(cells = self$cells) %>%
        left_join(celltype_map, by = c("cells", cell_colname)) %>%
        extract2(celltype_colname)
    },
    transcripts = function() {
      rownames(self$matrix)
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
      self$matrix %>%
        t() %>%
        as.matrix() %>%
        rowsum(group = self$celltypes) %>%
        t() %>%
        normalize_count_mat(
          type = self$params$normalization$type,
          scale = self$params$normalization$scale_factor
        )
    }
  )
)
