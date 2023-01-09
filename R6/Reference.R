library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Reference.R")
)

source(here("R6/CountMatrixWrapper.R"))

source(here("functions/norm_functions.R"))
source(here("R6/HeatmapPlot.R"))

Reference <- R6Class(
  "Reference",
  inherit = CountMatrixWrapper,
  public = list(
    initialize = function(markers, params) {
      private$.check_markers(markers)
      private$.check_params(params)

      super$initialize(markers$scrna_experiment$count_matrix)
      private$.params <- params
      private$.markers <- markers
    },
    print_heatmap = function(dir) {
      filename <- paste0(
        "thresh_", self$params$threshold,
        if (!is.null(self$params$cancer_celltypes)) {
          paste0(
            "-cancercols_",
            paste0(
              self$params$cancer_celltypes,
              collapse = "_"
            )
          )
        },
        ".png"
      )

      self$heatmap_plot$save(dir, filename)
    }
  ),
  active = list(
    count_matrix = function() {
      private$.count_matrix
    },
    markers = function() {
      private$.markers$threshold_quality(self$params$threshold)
    },
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
      if (is.null(private$.sigmat)) {
        private$.compute_sigmat()
      }
      private$.sigmat
    },
    heatmap_plot = function() {
      if (is.null(private$.heatmap_plot)) {
        private$.compute_heatmap_plot()
      }
      private$.heatmap_plot
    }
  ),
  private = list(
    .count_matrix = NULL,
    .markers = NULL,
    .params = NULL,
    .sigmat = NULL,
    .heatmap_plot = NULL,
    .compute_sigmat = function() {
      # TODO Check if this works as intended wrt normalization
      private$.sigmat <- self$count_matrix$celltype_count_matrix %>%
        magrittr::extract(
          i = rownames(.) %in% self$markers,
          j = !(colnames(.) %in% self$params$cancer_celltypes),
          drop = FALSE
        ) %>%
        normalize_count_mat(
          type = self$params$normalization$type,
          scale = self$params$normalization$scale_factor
        )
    },
    .compute_heatmap_plot = function() {
      title <- paste("Threshold:", self$params$threshold)

      private$.heatmap_plot <- HeatmapPlot$new(
        reference_matrix = self$sigmat,
        title = title
      )
    },
    .check_params = function(params) {
      stopifnot(
        all(c("R6", "Params", "ReferenceParams") %in% class(params))
      )
    },
    .check_markers = function(markers) {
      stopifnot(
        all(
          c("R6", "Markers") %in%
            class(markers)
        )
      )
    }
  )
)
