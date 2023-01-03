library("here")
library("R6")

suppressMessages(
  here::i_am("R6/DeconvolutionSummary.R")
)

source(here("functions/rmse.R"))

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
      # TODO Change to lazy field.
      self$list %>%
        lapply(
          function(deconv) {
            deconv$rmse
          }
        ) %>%
        unlist()
    },
    celltype_props_predicted_df = function() {
      if (is.null(private$.celltype_props_predicted_df)) {
        private$.compute_celltype_props_predicted_df()
      }
      private$.celltype_props_predicted_df
    },
    celltype_props_true_df = function() {
      if (is.null(private$.celltype_props_true_df)) {
        private$.compute_celltype_props_true_df()
      }
      private$.celltype_props_true_df
    },
    celltype_rmse_df = function() {
      if (is.null(private$.celltype_rmse_df)) {
        private$.compute_celltype_rmse_df()
      }
      private$.celltype_rmse_df
    }
  ),
  private = list(
    .params = NULL,
    .deconvolution_list = NULL,
    .param_df = NULL,
    .rmse_df = NULL,
    .rmse_summary = NULL,
    .celltype_props_predicted_df = NULL,
    .celltype_props_true_df = NULL,
    .celltype_rmse_df = NULL,
    .compute_param_df = function() {
      private$.param_df <- self$list %>%
        lapply(
          function(deconv) {
            deconv$params %>%
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
            merge(
              df1,
              df2,
              sort = FALSE,
              all = TRUE
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
    .compute_celltype_props_predicted_df = function() {
      celltype_props_predicted_df <- self$list %>%
        lapply(
          function(deconv) {
            deconv$deconvolution_output$proportions
          }
        ) %>%
        Reduce(
          f = function(df1, df2) {
            merge(
              df1,
              df2,
              sort = FALSE,
              all = TRUE
            )
          }
        )

      private$.celltype_props_predicted_df <- bind_cols(
        self$param_df,
        celltype_props_predicted_df
      )
    },
    .compute_celltype_props_true_df = function() {
      celltype_props_true_df <- self$list %>%
        lapply(
          function(deconv) {
            deconv$pseudobulk$celltype_abundances %>%
              t() %>%
              as.data.frame()
          }
        ) %>%
        Reduce(
          f = function(df1, df2) {
            merge(
              df1,
              df2,
              sort = FALSE,
              all = TRUE
            )
          }
        )

      private$.celltype_props_true_df <- bind_cols(
        self$param_df,
        celltype_props_true_df
      )
    },
    .compute_celltype_rmse_df = function() {
      # FIXME Provide a more sensible way to get celltypes.
      ctypes <- setdiff(
        names(self$celltype_props_true_df),
        self$param_vec
      )

      true_df <- self$celltype_props_true_df %>%
        bind_cols(data.frame(source = "true"))

      predicted_df <- self$celltype_props_predicted_df %>%
        bind_cols(data.frame(source = "predicted"))

      celltype_prop_df <- bind_rows(
        true_df,
        predicted_df
      ) %>%
        pivot_longer(
          all_of(ctypes),
          names_to = "celltype",
          values_to = "abundance"
        )

      if (any(duplicated(celltype_prop_df))) {
        warning(paste(
          "Detected & deduplicated non-unique parameter combinations while",
          "computing per celltype RMSEs. Check your deconvolution inputs if",
          "this is unexpected."
        ))
        celltype_prop_df %<>%
          filter(!duplicated(.))
      }

      private$.celltype_rmse_df <- celltype_prop_df %>%
        pivot_wider(
          id_cols = all_of(c(self$param_vec, "celltype")),
          names_from = source,
          values_from = abundance
        ) %>%
        # Removes NAs due to cancer_celltypes.
        drop_na() %>%
        group_by(across(all_of(c(self$params$summary_params, "celltype")))) %>%
        summarize(
          rmse = rmse(true, predicted)
        )
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
