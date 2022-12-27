library("here")
library("R6")

suppressMessages(
  here::i_am("R6/ReferenceParams.R")
)

source(here("R6/Params.R"))

ReferenceParams <- R6Class(
  "ReferenceParams",
  inherit = Params,
  public = list(
    metric = NULL,
    threshold = NULL,
    # TODO Change this from logical to char vec giving cols to split.
    split_cancer = NULL,
    cancer_celltypes = NULL,
    normalization = list(
      type = NULL,
      scale_factor = NULL
    ),
    initialize = function(metric = "DESeq2",
                          threshold = 0.5,
                          split_cancer = TRUE,
                          cancer_celltypes,
                          norm_type = "lognorm",
                          norm_scale_factor = 1) {
      param_list <- list(
        metric = metric,
        threshold = threshold,
        split_cancer = split_cancer,
        cancer_celltypes = cancer_celltypes,
        normalization = list(
          type = norm_type,
          scale_factor = norm_scale_factor
        )
      )

      private$.check_param_list(param_list)

      for (param_name in names(param_list)) {
        self[[param_name]] <- param_list[[param_name]]
      }
    }
  ),
  private = list(
    .check_param_list = function(param_list) {
      with(
        param_list,
        stopifnot(
          is.character(metric),
          length(metric) == 1,
          metric %in% c("deseq2"),
          is.numeric(threshold),
          length(threshold) == 1,
          threshold >= 0, threshold <= 1,
          is.logical(split_cancer),
          length(split_cancer) == 1,
          if (!is.null(cancer_celltypes)) is.character(cancer_celltypes),
          length(normalization) == 2,
          is.character(normalization$type),
          length(normalization$type) == 1,
          normalization$type %in% c("lognorm", "norm", "feature_scale"),
          is.numeric(normalization$scale_factor),
          length(normalization$scale_factor) == 1,
          normalization$scale_factor > 0
        )
      )
    }
  )
)
