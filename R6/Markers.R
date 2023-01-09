library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Markers.R")
)

# For the purpose of this project, Markers are a data.frame with (at least)
# two columns: transcript and metric. The former is just a character string
# identifying a transcript in the underlying ScRnaExperiment. The latter is
# some numeric measure that quantifies the fit of the transcript for being
# used in a reference.

Markers <- R6Class(
  "Markers",
  public = list(
    initialize = function(scrna_experiment) {
      private$.check_scrna_experiment(scrna_experiment)

      private$.scrna_experiment <- scrna_experiment
    },
    threshold_quality = function(qual_thresh) {
      stop(paste(
        "Method for thresholding by quality not implemented for Markers",
        "superclass."
      ))
    },
    threshold_quantity = function(quant_thresh = NULL, quant_prop = 0.1) {
      stop(paste(
        "Method for thresholding by quantity not implemented for Markers",
        "superclass."
      ))
    }
  ),
  active = list(
    scrna_experiment = function() {
      private$.scrna_experiment
    },
    matrix = function() {
      self$scrna_experiment$matrix_orig
    },
    meta = function() {
      self$scrna_experiment$meta
    },
    marker_df = function() {
      if (is.null(private$.marker_df)) {
        private$.compute_marker_df()
      }
      private$.marker_df
    }
  ),
  private = list(
    .scrna_experiment = NULL,
    .marker_df = NULL,
    .is_uniform = function(x) {
      # Test whether all elements of x are the same.
      all(x == x[1])
    },
    .compute_marker_df = function() {
      stop(paste(
        "Method for computing markers not implemented for Markers",
        "superclass."
      ))
    },
    .check_scrna_experiment = function(scrna_experiment) {
      stopifnot(
        all(
          c("R6", "CountMatrixWrapper", "ScRnaExperiment") %in%
            class(scrna_experiment)
        )
      )
    }
  )
)
