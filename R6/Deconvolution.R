library("deconvR")
library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Deconvolution.R")
)

soure(here("functions/rmse.R"))

Deconvolution <- R6Class(
  "Deconvolution",
  public = list(
    initialize = function(reference, pseudobulk, method) {
      stopifnot(
        all(c("Reference", "R6") %in% class(reference)),
        all(c("Pseudobulk", "R6") %in% class(pseudobulk))
      )
      private$.reference <- reference
      private$.pseudobulk <- pseudobulk
      private$.method <- method
    }
  ),
  active = list(
    reference = function() {
      private$.reference
    },
    pseudobulk = function() {
      private$.pseudobulk
    },
    method = function() {
      private$.method
    },
    res = function() {
      if (is.null(private$.celltype_props_predicted)) {
        private$.deconvolute()
      }
      private$.celltype_props_predicted
    },
    rmse = function() {
      if (is.null(private$.rmse)) {
        private$.compute_rmse()
      }
      privat$.rmse
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
    .reference = NULL,
    .pseudobulk = NULL,
    .method = NULL,
    .celltype_props_predicted = NULL,
    .rmse = NULL,
    .residuals_marker = NULL,
    .residuals_all = NULL,
    .deconvolute = function() {
      if (self$reference$n_celltypes <= 1) {
        # Some deconvolution methods will fail if only a single celltype is
        # present. This ensures the correct output.
        deconv_props <- 1 %>%
          set_names(self$reference$celltypes)
      } else {
        capture.output(
          suppressMessages(
            private$.celltype_props_predicted <- deconvR::deconvolute(
              # TODO Think about whether / how to ensure same transcripts in
              # ref & pbulk
              reference = self$reference$df,
              bulk = self$pseudobulk$df,
              model = self$method
            )$proportions
          ),
          type = c("output")
        )
      }
    },
    .compute_rmse = function() {
      rmse(
        private$.pseudobulk$celltype_abundances,
        self$.celltype_props_predicted
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
        self$.celltype_props_predicted

      transcript_abundances_marker <- pbulk$transcript_abundances %>%
        magrittr::extract(
          names(.) %in% ref$markers
        )

      private$.compute_residuals_marker <- transcript_abundances_marker -
        transcript_predictions_marker
    },
    .compute_residuals_all = function() {
      pbulk <- private$.pseudobulk
      ref <- private$.reference

      matrix_all <- pbulk$celltype_matrix_clean

      transcript_predictions_all <- matrix_all %*%
        self$.celltype_props_predicted

      transcript_abundances_all <- pbulk$transcript_abundances

      private$.compute_residuals_all <- transcript_abundances_allr -
        transcript_predictions_all
    }
  )
)
