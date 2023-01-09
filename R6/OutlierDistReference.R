library("here")
library("R6")

suppressMessages(
  here::i_am("R6/OutlierDistReference.R")
)

source(here("R6/Reference.R"))

OutlierDistReference <- R6Class(
  "OutlierDistReference",
  inherit = Reference,
  public = list(
    initialize = function(markers, params) {
      super$initialize(markers, params)
    }
  ),
  private = list(
    .compute_heatmap_plot = function() {
      title <- paste("Threshold:", self$params$threshold)

      private$.heatmap_plot <- HeatmapPlot$new(
        reference_matrix = self$sigmat,
        title = title,
        feature_scale = FALSE
      )
    },
    .check_markers = function(markers) {
      stopifnot(
        all(
          c("R6", "Markers", "OutlierDistMarkers") %in%
            class(markers)
        )
      )
    }
  )
)
