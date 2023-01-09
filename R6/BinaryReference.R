library("here")
library("R6")

here::i_am("R6/BinaryReference.R")

source(here("R6/Reference.R"))

BinaryReference <- R6Class(
  "BinaryReference",
  inherit = Reference,
  public = list(
    initialize = function(count_mat, threshold) {
      stop("Not implemented.")
    }
  )
)
