library("dplyr")
library("here")
library("R6")

here::i_am("R6/ScRnaExperiment.R")

source(here("R6/SigmatUtils"))
source(here("R6/CountMatrix"))
source(here("R6/Reference"))

ScRnaExperiment <- R6Class(
  "ScRnaExperiment",
  inherit = CountMatrix,
  public(
    initialize = function(count_mat_file, rowname_file, colname_file,
                          meta_file,
                          cell_col = "V1", celltype_col = "celltype_major",
                          downsample_frac = NULL) {
      private$sigmat_utils <- SigmatUtils$new()
      rename_vec <- c(
        "cell" = cell_col,
        "celltype" = celltype_col
      )

      meta <- fread(meta_file) %>%
        # See rlang::`!!!`.
        rename(!!!rename_vec)

      transcripts <- readLines(rowname_file)
      cells <- readLines(colname_file)

      matrix <- readMM(count_mat_file)

      matrix@Dimnames <- list(
        transcripts,
        cells
      )

      super$init(matrix, meta, downsample_frac)
    },
    create_reference = function(metric = c("raw_counts", "DESeq2"),
                                threshold) {
      # Metrics:
      #   "raw_counts": Produces a binary sigmat
      #   "DESeq2": *Selects* transcripts based on adjusted pvalue for
      #       differential expression, as determined by DESeq2.
      metric <- match.arg(metric)

      switch(metric,
        raw_counts = private$reference_raw_counts(threshold),
        DESeq2 = reference_DESeq2(threshold)
      )
    }
  ),
  private = list(
    sigmat_utils = NULL,
    # Marker dataframe is defined as a df with two cols:
    # metric and transcript.
    marker_df = NULL,
    get_marker_df = function() {
      if (is.null(private$marker_df)) {
        private$marker_df <-
          private$sigmatutils$create_marker_df(private$matrix)
      }
      private$marker_df
    },
    reference_raw_counts = function(threshold) {
      stop("Raw counts reference not implemented.")
    },
    reference_DESeq2 = function(threshold) {
      marker_thresh <- self$get_marker_df %>%
        slice_max(metric, prop = threshold)

      matrix <- private$matrix %>%
        extract(
          i = self$transcripts %in% marker_df_thresh$transcript,
          j = , drop = FALSE
        )

      meta <- private$meta %>%
        filter(rownames())
      Reference$new()
    }
  )
)
