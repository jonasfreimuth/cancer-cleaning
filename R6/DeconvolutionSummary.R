library("here")
library("R6")

suppressMessages(
  here::i_am("R6/DeconvolutionSummary.R")
)

DeconvolutionSummary <- R6Class(
  "DeconvolutionSummary",
  public = list(
    initialize = function(deconvolution_list) {
      private$.check_deconv_list(deconvolution_list)
      private$.deconvolution_list <- deconvolution_list
    }
  ),
  active = list(
    df = function() {
      stop("Not implemented.")
    },
    list = function() {
      private$.deconvolution_list
    }
  ),
  private = list(
    .deconvolution_list = NULL,
    .check_deconv_list = function(deconv_list) {
      all_deconvolution <- lapply(
        deconv_list,
        function(deconv) {
          all(c("Deconvolution", "R6") %in% class(deconv))
        }
      ) %>%
        unlist() %>%
        all()

      stopifnot(
        all_deconvolution
      )
    }
  )
)
