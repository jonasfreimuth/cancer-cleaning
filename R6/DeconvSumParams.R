library("here")
library("R6")

suppressMessages(
  here::i_am("R6/DeconvSumParams.R")
)

source(here("R6/Params.R"))

DeconvSumParams <- R6Class(
  "DeconvSumParams",
  inherit = Params,
  public = list(
    .param_list = NULL,
    summary_params = NULL,
    initialize = function(summary_params = c(
                            "reference.threshold", "reference.cancer_celltypes"
                          )) {
      param_list <- list(
        summary_params = summary_params
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
          is.character(summary_params)
        )
      )
    }
  )
)
