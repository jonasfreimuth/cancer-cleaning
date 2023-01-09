library("dplyr")
library("here")
library("R6")

suppressMessages(
  here::i_am("R6/ScRnaExperiment.R")
)

source(here("R6/SigmatUtils.R"))
source(here("R6/CountMatrix.R"))
source(here("R6/CountMatrixWrapper.R"))
source(here("R6/Deseq2Markers.R"))
source(here("R6/Deseq2Reference.R"))
source(here("R6/OutlierDistMarkers.R"))
source(here("R6/OutlierDistReference.R"))
source(here("R6/Pseudobulk.R"))

ScRnaExperiment <- R6Class(
  "ScRnaExperiment",
  inherit = CountMatrixWrapper,
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

      super$initialize(CountMatrix$new(matrix, meta, downsample_frac))
    },
    create_reference = function(params) {
      stopifnot(
        all(c("R6", "Params", "ReferenceParams") %in% class(params))
      )
      # Metrics:
      #   "raw_counts": Produces a binary sigmat
      #   "deseq2": *Selects* transcripts based on adjusted pvalue for
      #       differential expression, as determined by DESeq2.
      switch(params$metric,
        raw_counts = private$.reference_raw_counts(params),
        deseq2 = private$.reference_deseq2(params),
        outlier_dist = private$.reference_outlier_dist(params)
      )
    },
    create_pseudobulk = function(params) {
      stopifnot(
        all(c("R6", "Params", "PseudobulkParams") %in% class(params))
      )
      Pseudobulk$new(private$.count_matrix, params)
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
    .marker_df = NULL,
    .get_marker_df = function() {
      if (is.null(private$.marker_df)) {
        private$.marker_df <-
          private$.sigmat_utils$create_marker_df(
            self$matrix_orig,
            self$meta
          )
      }
      private$.marker_df
    },
    .reference_raw_counts = function(params) {
      stop("Raw counts reference not implemented.")
    },
    .reference_deseq2 = function(params) {
      Deseq2Reference$new(
        Deseq2Markers$new(
          self
        ),
        params
      )
    },
    .reference_outlier_dist = function(params) {
      OutlierDistReference$new(
        OutlierDistMarkers$new(
          self
        ),
        params
      )
    }
  )
)
