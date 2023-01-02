library("here")
library("R6")

suppressMessages(
  here::i_am("R6/PseudobulkParams.R")
)

source(here("R6/Params.R"))

PseudobulkParams <- R6Class(
  "PseudobulkParams",
  inherit = Params,
  public = list(
    .param_list = NULL,
    id = NULL,
    cell_indices = NULL,
    cancer_celltypes = NULL,
    normalization = list(
      type = NULL,
      scale_factor = 1
    ),
    initialize = function(id,
                          cell_indices,
                          cancer_celltypes = "",
                          norm_type = NULL,
                          norm_scale_factor = 1) {
      if (is.null(cancer_celltypes)) {
        cancer_celltypes <- ""
      }

      param_list <- list(
        id = as.integer(id),
        cell_indices = as.integer(cell_indices),
        cancer_celltypes = cancer_celltypes,
        normalization = list(
          type = norm_type,
          scale_factor = norm_scale_factor
        )
      )

      private$.check_param_list(param_list)

      self$.param_list <- param_list

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
          is.integer(id),
          length(id) == 1,
          id >= 0,
          is.integer(cell_indices),
          length(cell_indices) > 0,
          !any(duplicated(cell_indices)),
          is.character(cancer_celltypes),
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
