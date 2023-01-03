library("deconvR")
library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Deconvolution.R")
)

source(here("R6/DeconvParams.R"))

source(here("functions/rmse.R"))

Deconvolution <- R6Class(
  "Deconvolution",
  public = list(
    initialize = function(reference, pseudobulk, params) {
      private$.check_params(params)
      private$.params <- params

      stopifnot(
        all(c("Reference", "R6") %in% class(reference)),
        all(c("Pseudobulk", "R6") %in% class(pseudobulk))
      )
      private$.check_ref_pbulk_fit(reference, pseudobulk)
      private$.reference <- reference
      private$.pseudobulk <- pseudobulk
    }
  ),
  active = list(
    # TODO Handle params differently so that it fits the general template (i.e.
    # params active returns just the params of the object.)
    deconv_params = function() {
      private$.params
    },
    params = function() {
      list(
        reference = self$reference$params$.param_list,
        pseudobulk = self$pseudobulk$params$.param_list,
        deconvolution = self$deconv_params$.param_list
      )
    },
    reference = function() {
      private$.reference
    },
    pseudobulk = function() {
      private$.pseudobulk
    },
    method = function() {
      private$.method
    },
    deconvolution_output = function() {
      if (is.null(private$.deconvolution_output)) {
        private$.compute_deconvolution_output()
      }
      private$.deconvolution_output
    },
    celltype_props_predicted = function() {
      self$deconvolution_output$proportions %>%
        as.matrix()
    },
    rmse = function() {
      if (is.null(private$.rmse)) {
        private$.compute_rmse()
      }
      private$.rmse
    },
    residuals_marker = function() {
      if (is.null(private$.residuals_marker)) {
        private$.compute_residuals_marker()
      }
      private$.residuals_marker
    },
    residuals_all = function() {
      if (is.null(private$.residuals_all)) {
        private$.compute_residuals_all()
      }
      private$.residuals_all
    }
  ),
  private = list(
    .params = NULL,
    .reference = NULL,
    .pseudobulk = NULL,
    .method = NULL,
    .deconvolution_output = NULL,
    .rmse = NULL,
    .residuals_marker = NULL,
    .residuals_all = NULL,
    .cor_marker = NULL,
    .cor_all = NULL,
    .compute_deconvolution_output = function() {
      if (self$reference$n_celltypes <= 1) {
        # Some deconvolution methods will fail if only a single celltype is
        # present. This ensures the correct output.
        deconv_props <- 1 %>%
          set_names(self$reference$celltypes)
      } else {
        capture.output(
          suppressMessages(
            private$.deconvolution_output <- deconvR::deconvolute(
              # TODO Think about whether / how to ensure same transcripts in
              # ref & pbulk
              reference = self$reference$df,
              bulk = self$pseudobulk$df,
              model = self$deconv_params$deconvolution_method
            )
          ),
          type = c("output")
        )
      }
    },
    .compute_rmse = function() {
      private$.rmse <- rmse(
        self$pseudobulk$celltype_abundances,
        self$celltype_props_predicted
      )
    },
    .compute_residuals_marker = function() {
      pbulk <- private$.pseudobulk
      ref <- private$.reference

      matrix_marker <- pbulk$celltype_matrix_clean %>%
        magrittr::extract(
          i = rownames(.) %in% ref$markers,
          j = ,
          drop = FALSE
        )

      transcript_predictions_marker <- matrix_marker %*%
        t(self$celltype_props_predicted)

      transcript_abundances_marker <- pbulk$transcript_abundances %>%
        magrittr::extract(
          names(.) %in% ref$markers
        )

      private$.residuals_marker <- transcript_abundances_marker -
        transcript_predictions_marker
    },
    .compute_residuals_all = function() {
      pbulk <- private$.pseudobulk
      ref <- private$.reference

      matrix_all <- pbulk$celltype_matrix_clean

      transcript_predictions_all <- matrix_all %*%
        t(self$celltype_props_predicted)

      transcript_abundances_all <- pbulk$transcript_abundances

      private$.residuals_all <- transcript_abundances_all -
        transcript_predictions_all
    },
    .compute_cor_marker = function() {
      # TODO Trim transcript_abundances_cancer to residuals
      private$.cor_marker <- cor(
        self$residuals_marker,
        private$.pseudobulk$transcript_abundances_cancer,
        method = self$deconv_params$correlation_method
      )
    },
    .compute_cor_all = function() {
      private$.cor_all <- cor(
        self$residuals_all,
        private$.pseudobulk$transcript_abundances_cancer,
        method = self$deconv_params$correlation_method
      )
    },
    .check_params = function(params) {
      stopifnot(
        all(c("R6", "Params", "DeconvParams") %in% class(params))
      )
    },
    .check_ref_pbulk_fit = function(reference, pseudobulk) {
      ref_cancer_ctypes <- reference$params$cancer_celltypes
      pbulk_cancer_ctypes <- pseudobulk$params$cancer_celltypes

      stopifnot(
        isTRUE(all.equal(sort(ref_cancer_ctypes), sort(pbulk_cancer_ctypes)))
      )
    }
  )
)
