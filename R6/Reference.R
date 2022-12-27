library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Reference.R")
)

source(here("R6/CountMatrixWrapper.R"))

source(here("functions/norm_functions.R"))

Reference <- R6Class(
  "Reference",
  inherit = CountMatrixWrapper,
  public = list(
    initialize = function(count_matrix, markers, params) {
      # TODO Check more args.
      private$.check_markers(markers)
      private$.check_count_matrix(count_matrix)
      private$.check_params(params)

      super$initialize(count_matrix)
      private$.params <- params
      private$.markers <- markers
    }
  ),
  active = list(
    params = function() {
      private$.params
    },
    matrix_raw = function() {
      # TODO Change to matrix_orig where matrix_raw is used.
      self$matrix_orig
    },
    meta_raw = function() {
      private$.count_matrix$meta
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
          i = rownames(.) %in% private$.markers,
          j = !(colnames(.) %in% self$params$cancer_celltypes),
          drop = FALSE
        ) %>%
        normalize_count_mat(
          type = self$params$nomalization$type,
          scale = self$params$nomalization$scale_factor
        )
    }
  ),
  private = list(
    .count_matrix = NULL,
    .markers = NULL,
    .params = NULL,
    .check_markers = function(markers) {
      stopifnot(
        is.character(markers)
      )
    },
    .check_params = function(params) {
      stopifnot(
        all(c("R6", "Params", "ReferenceParams") %in% class(params))
      )
    },
    .check_count_matrix = function(count_matrix) {
      stopifnot(
        all(c("R6", "CountMatrix") %in% class(count_matrix))
      )
    }
  )
)
