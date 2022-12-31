library("here")
library("R6")

here::i_am("R6/ReferenceBinary.R")

source(here("R6/Reference.R"))

ReferenceBinary <- R6Class(
  "ReferenceBinary",
  inherit = Reference,
  public = list(
    initialize = function(count_mat, threshold) {
      stop("Not implemented.")
    }
  )
)
