library("DESeq2")
library("scuttle")
library("dplyr")
library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Deseq2Markers.R")
)

source(here("R6/Markers.R"))

Deseq2Markers <- R6Class(
  "Deseq2Markers",
  inherit = Markers,
  public = list(
    initialize = function(scrna_experiment) {
      super$initialize(scrna_experiment)

      private$.design <- ~celltype
    },
    threshold_quality = function(qual_thresh) {
      # Thresholding based on pvalue of DE analysis for log fold change.
      self$marker_df %>%
        filter(metric <= qual_thresh) %>%
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
        arrange(metric) %>%
        magrittr::extract(i = seq_len(quant_thresh))
    }
  ),
  active = list(
    reduced_design = function() {
      if (length(self$design) >= 3) {
        # This assumes the last term of the design formula is the only one of
        # interest.
        # TODO Enforce this assumption.
        return(self$design[-length(self$design)])
      } else {
        return(~1)
      }
    },
    design = function() {
      private$.design
    },
    size_factors = function() {
      private$.get_size_factors()
    },
    matrix_prefiltered = function() {
      if (is.null(private$.matrix_prefiltered)) {
        private$.compute_matrix_prefiltered()
      }
      private$.matrix_prefiltered
    },
    matrix_clean = function() {
      if (is.null(private$.matrix_clean)) {
        private$.compute_matrix_clean()
      }
      private$.matrix_clean
    },
    meta_clean = function() {
      if (is.null(private$.meta_clean)) {
        private$.compute_meta_clean()
      }
      private$.meta_clean
    },
    ds2_data = function() {
      if (is.null(private$.ds2_data)) {
        private$.compute_ds2_data()
      }
      private$.ds2_data
    },
    de_transcripts = function() {
      if (is.null(private$.de_transcripts)) {
        private$.compute_de_transcripts()
      }
      private$.de_transcripts
    },
    marker_df = function() {
      if (is.null(private$.marker_df)) {
        private$.compute_marker_df()
      }
      private$.marker_df
    }
  ),
  private = list(
    .matrix_clean = NULL,
    .meta_clean = NULL,
    .matrix_prefiltered = NULL,
    .design = NULL,
    .ds2_data = NULL,
    .de_transcripts = NULL,
    .compute_matrix_prefiltered = function() {
      private$.matrix_prefiltered <- self$matrix %>%
        # remove rows that are all the same
        # Also includes rows that are all zero-counts
        magrittr::extract(
          i = !apply(
            ., 1, private$.is_uniform
          ),
          j = , drop = FALSE
        ) %>%
        magrittr::extract(
          i = ,
          j = order(colnames(.)),
          drop = FALSE
        )
    },
    # TODO Check if this can be merged with matrix_prefiltered
    .compute_matrix_clean = function() {
      private$.matrix_clean <- self$matrix_prefiltered %>%
        magrittr::extract(
          i = ,
          j = !colSums(.) == 0,
          drop = FALSE
        )

      zero_cols <- setdiff(self$matrix_clean, self$matrix_prefiltered)
      if (length(zero_cols) >= 1) {
        warning(paste(
          "During calculation of size factors, cells",
          paste(zero_cols, collapse = ", "),
          "had no transcript",
          "expression after removal of uniform rows and were therefore not",
          "considered."
        ))
      }
    },
    .compute_meta_clean = function() {
      private$.meta_clean <- self$meta %>%
        # Prevent DESeq fun from complaining
        mutate(across(where(is.character), as.factor)) %>%
        filter(cell %in% colnames(self$matrix_clean))
    },
    .get_size_factors = function() {
      # TODO Is prescaling useful here?
      withCallingHandlers(
        self$matrix_prefiltered %>%
          scuttle::pooledSizeFactors(
            # TODO Think about whether there are general improvements possible
            # that would eliminate the necessity for positive = TRUE
            positive = TRUE
          ),
        warning = function(w) {
          is_size_factor_warning <- {
            w$message == "encountered non-positive size factor estimates"
          }
          if (is_size_factor_warning) {
            tryInvokeRestart("muffleWarning")
          }
        }
      )
    },
    .compute_ds2_data = function() {
      private$.ds2_data <- DESeqDataSetFromMatrix(
        countData = self$matrix_clean,
        colData = self$meta_clean,
        design = self$design
      ) %>%
        `sizeFactors<-`(value = self$size_factors)
    },
    .compute_de_transcripts = function() {
      private$.de_transcripts <- self$ds2_data %>%
        DESeq(
          # Args are set according to
          # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis.
          test = "LRT",
          fitType = "glmGamPoi",

          # This is recommended, but according to the source above shouldn't
          # do anything with test = "LRT".
          useT = TRUE,

          # This should also be already set when using fitType = "glmGamPoi".
          minmu = 10^-6,
          reduced = self$reduced_design,
          quiet = TRUE
        ) %>%
        # Caution, results may contain NAs, seems to mainly depend on whether
        # zero cols were cleaned beforehand.
        results(
          tidy = TRUE
        )
    },
    .compute_marker_df = function() {
      private$.marker_df <- self$de_transcripts %>%
        # This may be unncecessary if uniform rows have been previously
        # removed.
        drop_na(!lfcSE) %>%
        arrange(padj) %>%
        # Transform to generic form of a df with a transcript col and a metric
        # col, with the latter being some metric that can be thresholded to
        # select most informative marker transcripts.
        select(
          transcript = row,
          metric = padj
        )
    }
  )
)
