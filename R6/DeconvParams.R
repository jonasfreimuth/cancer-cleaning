library("here")
library("R6")

suppressMessages(
  here::i_am("R6/DeconvParams.R")
)

source(here("R6/Params.R"))

DeconvParams <- R6Class(
  "DeconvParams",
  inherit = Params,
  public = list(
    .param_list = NULL,
    deconvolution_method = NULL,
    correlation_method = NULL,
    initialize = function(deconvolution_method = "nnls",
                          correlation_method = "pearson") {
      param_list <- list(
        deconvolution_method = deconvolution_method,
        correlation_method = correlation_method
      )

      private$.check_param_list(param_list)

      self$.param_list <- param_list

      for (param_name in names(param_list)) {
        self[[param_name]] <- param_list[[param_name]]
      }
    },
    check_summary_params_fit = function(params) {
      stopifnot(
        all(self$summary_params %in% params)
      )
    }
  ),
  private = list(
    .check_param_list = function(param_list) {
      with(
        param_list,
        stopifnot(
          is.character(deconvolution_method),
          deconvolution_method %in% c("nnls", "qp", "rlm", "svm"),
          is.character(correlation_method),
          correlation_method %in% c("pearson", "kendall", "spearman")
        )
      )
    }
  )
)
