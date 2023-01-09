library("here")
library("R6")

suppressMessages(
  here::i_am("R6/OutlierDistMarkers.R")
)

source(here("R6/Markers.R"))

OutlierDistMarkers <- R6Class(
  "OutlierDistMarkers",
  inherit = Markers,
  public = list(
    initialize = function(scrna_experiment) {
      super$initialize(scrna_experiment)
    },
    threshold_quality = function(qual_thresh) {
      # Thresholding based on distance of outlier from the average distance.
      self$marker_df %>%
        filter(metric >= qual_thresh) %>%
        magrittr::extract2("transcript")
    },
    threshold_quantity = function(quant_thresh = NULL, quant_prop = 0.1) {
      stopifnot(
        xor(!is.null(quant_thresh), !is.null(quant_prop))
      )

      if (is.null(quant_thresh)) {
        quant_thresh <- quant_prop * nrow(self$marker_df)
      }

      self$marker_df %>%
        arrange(desc(metric)) %>%
        magrittr::extract(i = seq_len(quant_thresh))
    }
  ),
  active = list(
    matrix_fs = function() {
      if (is.null(private$.matrix_fs)) {
        private$.compute_matrix_fs()
      }
      private$.matrix_fs
    }
  ),
  private = list(
    .matrix_fs = NULL,
    .hampel_interval = function(x, mad_mult = 3) {
      x_median <- median(x)
      x_mad_mult <- mad(x) * mad_mult

      return(c(x_median - x_mad_mult, x_median + x_mad_mult))
    },
    .outlier_dist = function(x, mad_mult = 3) {
      dists <- mean_dists(x)
      dist_interval <- private$.hampel_interval(dists, mad_mult)

      outlier_idcs <- dists < dist_interval[1] | dists > dist_interval[2]

      if (sum(outlier_idcs) == 1) {
        return(dists[outlier_idcs])
      } else {
        return(0)
      }
    },
    .compute_matrix_fs = function() {
      private$.matrix_fs <- self$matrix %>%
        magrittr::extract(
          i = !apply(., 1, private$.is_uniform),
          j = , drop = FALSE
        ) %>%
        apply(
          1,
          feature_scale
        ) %>%
        t()
    },
    .compute_marker_df = function() {
      private$.marker_df <- self$matrix_fs %>%
        apply(
          1,
          private$.outlier_dist
        ) %>%
        {
          data.frame(
            transcript = names(.),
            metric = .
          )
        }
    }
  )
)
