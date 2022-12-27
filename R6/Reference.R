library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Reference.R")
)

Reference <- R6Class(
  "Reference",
  public = list(
    # TODO Consider getting params into private.
    params = list(
      metric = NULL,
      threshold = NULL,
      # TODO Change this from logical to char vec giving cols to split.
      split_cancer = NULL,
      cancer_celltypes = NULL,
      normalization = list(
        type = NULL,
        scale_factor = NULL
      )
    ),
    initialize = function(count_matrix, threshold, markers) {
      # TODO Check more args.
      private$.check_markers(markers)
      private$.check_threshold(threshold)
      private$.check_count_matrix(count_matrix)

      private$.count_matrix <- count_matrix
      self$params$threshold <- threshold
      private$.markers <- markers
    }
  ),
  active = list(
    matrix_raw = function() {
      private$.count_marix$matrix
    },
    meta_raw = function() {
      private$.count_matrix$meta
    },
    transcripts = function() {
      rownames(self$sigmat)
    },
    n_transcripts = function() {
      nrow(self$sigmat)
    },
    celltypes = function() {
      colnames(self$sigmat)
    },
    n_celltypes = function() {
      ncol(self$sigmat)
    },
    markers = function() {
      private$.markers
    },
    df = function() {
      self$sigmat %>%
        {
          # TODO Check whether this is necessary (sigmat might be sparse)
          as.matrix(.) %>%
            as.data.frame() %>%
            mutate(IDs = rownames(.))
        } %>%
        select(IDs, everything())
    },
    sigmat = function() {
      # TODO Change this to lazy field.
      # TODO Check if this works as intended wrt normalization
      private$.count_matrix$celltype_count_matrix %>%
        magrittr::extract(
          i = private$.markers,
          j = ,
          drop = FALSE
        )
    }
  ),
  private = list(
    .count_matrix = NULL,
    .markers = NULL,
    .check_markers = function(markers) {
      stopifnot(
        is.character(markers)
      )
    },
    .check_threshold = function(threshold) {
      stopifnot(
        is.numeric(threshold),
        length(threshold) == 1
      )
    },
    .check_count_matrix = function(count_matrix) {
      stopifnot(
        c("R6", "CountMatrix") %in% class(count_matrix)
      )
    }
  )
)
