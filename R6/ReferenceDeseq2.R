library("here")
library("R6")

here::i_am("R6/ReferenceDeseq2.R")

source(here("R6/Reference.R"))

ReferenceDeseq2 <- R6Class(
  "ReferenceDeseq2",
  inherit = Reference,
  public = list(
    initialize = function(count_mat, markers, params) {
      super$initialize(count_mat, markers, params)
    }
  )
)
