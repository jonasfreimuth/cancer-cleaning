library("deconvR")
library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Deconvolution.R")
)

source(here("functions/rmse.R"))

Deconvolution <- R6Class(
  "Deconvolution",
  public = list(
    initialize = function(reference, pseudobulk, method) {
      stopifnot(
        all(c("Reference", "R6") %in% class(reference)),
        all(c("Pseudobulk", "R6") %in% class(pseudobulk))
      )
      private$.check_ref_pbulk_fit(reference, pseudobulk)
      private$.reference <- reference
      private$.pseudobulk <- pseudobulk
      # TODO Turn this into a proper param object.
      private$.method <- method
    }
  ),
  active = list(
    params = function() {
      list(
        reference = self$reference$params$.param_list,
        pseudobulk = self$pseudobulk$params$.param_list,
        deconvolution = list(
          method = self$method
        )
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
    .reference = NULL,
    .pseudobulk = NULL,
    .method = NULL,
    .celltype_props_predicted = NULL,
    .rmse = NULL,
    .residuals_marker = NULL,
    .residuals_all = NULL,
    .cor_marker = NULL,
    .cor_all = NULL,
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
      private$.rmse <- rmse(
        self$pseudobulk$celltype_abundances,
        self$res
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
    },
    .compute_cor_marker = function() {
      method <- "pearson"
      # TODO Trim transcript_abundances_cancer to residuals
      private$.cor_marker <- cor(
        self$residuals_marker,
        private$.pseudobulk$transcript_abundances_cancer,
        method = method
      )
    },
    .compute_cor_all = function() {
      method <- "pearson"
      private$.cor_all <- cor(
        self$residuals_all,
        private$.pseudobulk$transcript_abundances_cancer,
        method = method
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
