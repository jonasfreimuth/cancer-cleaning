library("here")
library("R6")

here::i_am("R6/Deseq2Reference.R")

source(here("R6/Reference.R"))

Deseq2Reference <- R6Class(
  "Deseq2Reference",
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
          c("R6", "Markers", "Deseq2Markers") %in%
            class(markers)
        )
      )
    }
  )
)
