library("fastDummies")
library("DESeq2")
library("scuttle")
library("dplyr")
library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Deseq2TtestMarkers.R")
)

source(here("R6/Markers.R"))

Deseq2TtestMarkers <- R6Class(
  "Deseq2TtestMarkers",
  inherit = Markers,
  public = list(
    initialize = function(scrna_experiment) {
      super$initialize(scrna_experiment)

      private$.design <- ~celltype_dummy
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
      # For Deseq2Ttest, size factors are needed multiple times, so it's a lazy
      # field (as opposed to the simple active one in Deseq2).
      private$.compute_size_factors()
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
    dummy_df = function() {
      if (is.null(private$.dummy_df)) {
        private$.compute_dummy_df()
      }
      private$.dummy_df
    },
    ds2_data_list = function() {
      if (is.null(private$.ds2_data_list)) {
        private$.compute_ds2_data_list()
      }
      private$.ds2_data_list
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
    .dummy_df = NULL,
    .matrix_prefiltered = NULL,
    .design = NULL,
    .size_factors = NULL,
    .ds2_data_list = NULL,
    .de_transcripts = NULL,
    .compute_matrix_prefiltered = function() {
      private$.matrix_prefiltered <- self$matrix %>%
        # Prevent DESeq fun from complaining
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
    .compute_dummy_df = function() {
      private$.dummy_df <- self$meta_clean %>%
        select(celltype) %>%
        fastDummies::dummy_cols(
          select_columns = "celltype",
          remove_selected_columns = TRUE
        ) %>%
        rename_with(
          .fn = function(name) {
            str_replace_all(
              name,
              paste0("celltype", "_"),
              ""
            )
          }
        )
    },
    .compute_size_factors = function() {
      # TODO Is prescaling useful here?
      private$.size_factors <- withCallingHandlers(
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
    .compute_ds2_data_list = function() {
      # This is to just shut DESeqDataSetFromMatrix up. matrix_clean is still
      # sparse and if we can avoid it we shouldn't save non-sparse matrices
      # because they may be very big. However I couldn't find an easier way to
      # coerce it to integer as the DESeq... fun wants its input to be.
      matrix_int <- self$matrix_clean %>%
        as.matrix() %>%
        `mode<-`("integer")

      private$.ds2_data_list <- names(self$dummy_df) %>%
        set_names(.) %>%
        lapply(
          function(celltype_col) {
            dummy_meta <- self$meta_clean %>%
              select(cell) %>%
              mutate(
                celltype_dummy = as.factor(self$dummy_df[[celltype_col]])
              )

            ds2_data <- DESeqDataSetFromMatrix(
              countData = matrix_int,
              colData = dummy_meta,
              design = self$design
            ) %>%
              `sizeFactors<-`(value = self$size_factors)

            return(ds2_data)
          }
        )
    },
    .compute_de_transcripts = function() {
      private$.de_transcripts <- self$ds2_data_list %>%
        lapply(
          function(ds2_data) {
            DESeq(
              ds2_data,
              # Args are set according to
              # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis.
              test = "LRT",
              fitType = "glmGamPoi",

              # This is recommended, but according to the source above shouldn't
              # do anything with test = "LRT".
              useT = TRUE,

              # This should also be already set when using
              # fitType = "glmGamPoi".
              minmu = 10^-6,
              reduced = self$reduced_design,
              quiet = TRUE
            ) %>%
              # Caution, results may contain NAs, seems to mainly depend on
              # whether zero cols were cleaned beforehand.
              results(
                tidy = TRUE
              )
          }
        ) %>%
        bind_rows(.id = "celltype")
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
