library("dplyr")
library("here")
library("R6")

here::i_am("R6/ScRnaExperiment.R")

source(here("R6/SigmatUtils.R"))
source(here("R6/CountMatrix.R"))
source(here("R6/Reference.R"))

ScRnaExperiment <- R6Class(
  "ScRnaExperiment",
  public = list(
    initialize = function(count_mat_file, rowname_file, colname_file,
                          meta_file,
                          cell_col = "V1", celltype_col = "celltype_major",
                          downsample_frac = NULL) {
      private$.sigmat_utils <- SigmatUtils$new()

      .rename_fun <- function(x) {
        recode_vec <- c("cell", "celltype") %>%
          set_names(cell_col, celltype_col)

        # See rlang::`!!!`.
        dplyr::recode(x, !!!recode_vec)
      }

      meta <- fread(meta_file) %>%
        dplyr::rename_with(.fn = .rename_fun)

      transcripts <- readLines(rowname_file)
      cells <- readLines(colname_file)

      matrix <- readMM(count_mat_file)

      matrix@Dimnames <- list(
        transcripts,
        cells
      )

      private$.count_matrix <- CountMatrix$new(matrix, meta, downsample_frac)
    },
    create_reference = function(metric = c("raw_counts", "DESeq2"),
                                threshold) {
      # Metrics:
      #   "raw_counts": Produces a binary sigmat
      #   "DESeq2": *Selects* transcripts based on adjusted pvalue for
      #       differential expression, as determined by DESeq2.
      metric <- match.arg(metric)

      switch(metric,
        raw_counts = private$.reference_raw_counts(threshold),
        DESeq2 = reference_DESeq2(threshold)
      )
    }
  ),
  active = list(
    count_matrix = function() {
      private$.count_matrix
    }
  ),
  private = list(
    .count_matrix = NULL,
    .sigmat_utils = NULL,
    # Marker dataframe is defined as a df with two cols:
    # metric and transcript.
    .marker_df = NULL,
    .get_marker_df = function() {
      if (is.null(private$.marker_df)) {
        private$.marker_df <-
          private$.sigmatutils$create_marker_df(private$.count_matrix$matrix)
      }
      private$.marker_df
    },
    .reference_raw_counts = function(threshold) {
      stop("Raw counts reference not implemented.")
    },
    .reference_DESeq2 = function(threshold) {
      marker_thresh <- self$get_marker_df %>%
        slice_max(metric, prop = threshold)

      reference_transcripts <- marker_thresh$transcript

      Reference$new(private$.count_matrix, reference_transcripts, threshold)
    }
  )
)
