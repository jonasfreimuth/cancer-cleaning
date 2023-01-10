library("here")
library("R6")

here::i_am("R6/Deseq2TtestReference.R")

source(here("R6/Reference.R"))

Deseq2TtestReference <- R6Class(
  "Deseq2TtestReference",
  inherit = Reference,
  public = list(
    initialize = function(markers, params) {
      super$initialize(markers, params)
    }
  ),
  private = list(
    .check_markers = function(markers) {
      stopifnot(
        all(
          c("R6", "Markers", "Deseq2TtestMarkers") %in%
            class(markers)
        )
      )
    }
  )
)
