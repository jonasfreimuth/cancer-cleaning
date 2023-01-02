library("here")
library("R6")

suppressMessages(
  here::i_am("R6/DeconvolutionSummary.R")
)

DeconvolutionSummary <- R6Class(
  "DeconvolutionSummary",
  public = list(
    initialize = function(deconvolution_list,
                          params) {
      private$.check_deconv_list(deconvolution_list)
      private$.check_params(params)

      private$.deconvolution_list <- deconvolution_list
      private$.params <- params

      self$params$check_summary_params_fit(self$param_vec)
    }
  ),
  active = list(
    params = function() {
      private$.params
    },
    param_df = function() {
      if (is.null(private$.param_df)) {
        private$.compute_param_df()
      }
      private$.param_df
    },
    param_vec = function() {
      names(self$param_df)
    },
    df = function() {
      stop("Not implemented.")
    },
    list = function() {
      private$.deconvolution_list
    },
    rmse_summary = function() {
      if (is.null(private$.rmse_summary)) {
        private$.compute_rmse_summary()
      }
      private$.rmse_summary
    },
    rmse_df = function() {
      if (is.null(private$.rmse_df)) {
        private$.compute_rmse_df()
      }
      private$.rmse_df
    },
    rmse_vec = function() {
      self$list %>%
        lapply(
          function(deconv) {
            deconv$rmse
          }
        ) %>%
        unlist()
    }
  ),
  private = list(
    .params = NULL,
    .deconvolution_list = NULL,
    .param_df = NULL,
    .rmse_df = NULL,
    .rmse_summary = NULL,
    .compute_param_df = function() {
      private$.param_df <- self$list %>%
        lapply(
          function(deconv) {
            ref_params <- deconv$params$reference
            pbulk_params <- deconv$params$pseudobulk

            list(
              reference = ref_params,
              pseudobulk = pbulk_params
            ) %>%
              rapply(
                paste,
                classes = "ANY",
                how = "replace",
                collapse = ";"
              ) %>%
              unlist() %>%
              t() %>%
              as.data.frame()
          }
        ) %>%
        Reduce(
          f = function(df1, df2) {
            merge(df1,
              df2,
              by = intersect(colnames(df1), colnames(df2)),
              all.x = TRUE,
              all.y = TRUE
            )
          }
        )
    },
    .compute_rmse_df = function() {
      rmse_df <- self$param_df
      rmse_df$rmse <- self$rmse_vec
      private$.rmse_df <- rmse_df
    },
    .compute_rmse_summary = function() {
      private$.rmse_summary <- self$rmse_df %>%
        group_by(across(self$params$summary_params)) %>%
        summarize(mean_rmse = mean(rmse), .groups = "drop_last")
    },
    .check_params = function(params) {
      stopifnot(
        all(c("R6", "Params", "DeconvSumParams") %in% class(params))
      )
    },
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
