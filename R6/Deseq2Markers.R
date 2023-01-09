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
  private = list(
    .compute_marker_df = function() {
      matrix <- self$scrna_experiment$matrix_orig
      meta <- self$scrna_experiment$meta

      private$.marker_df <- matrix %>%
        # remove rows that are all the same
        # Also includes rows that are all zero-counts
        magrittr::extract(
          i = !apply(
            ., 1, private$.is_uniform
          ), j = , drop = FALSE
        ) %>%
        private$.get_de_transcripts(
          meta,
          ~celltype
        ) %>%
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
    },
    .get_de_transcripts = function(count_mat, meta, design) {
      # TODO Get the functions in here to STFU.
      meta <- meta %>%
        arrange(cell) %>%
        # Prevent DESeq fun from complaining
        mutate(across(where(is.character), as.factor))
      count_mat <- count_mat[, order(colnames(count_mat))]

      zero_sum_cells <- colSums(count_mat) == 0

      if (any(zero_sum_cells)) {
        # TODO Make this a more helpfull message.
        warning(paste(
          "During calculation of size factors, cells",
          paste(colnames(count_mat)[zero_sum_cells], collapse = ", "),
          "had no transcript",
          "expression after removal of uniform rows and were not considered."
        ))
        count_mat <- count_mat[, !zero_sum_cells]
        meta <- meta[!zero_sum_cells, ]
      }

      # TODO Is prescaling useful here?
      size_factors <- count_mat %>%
        scuttle::pooledSizeFactors()

      ds2_data <- DESeqDataSetFromMatrix(
        countData = count_mat,
        colData = meta,
        design = design
      )

      if (length(design) >= 3) {
        # This assumes the last term of the design formula is the only one of
        # interest.
        # TODO Enforce this assumption.
        reduced_design <- design[-length(design)]
      } else {
        reduced_design <- ~1
      }

      ds2_data %>%
        `sizeFactors<-`(value = size_factors) %>%
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
          reduced = reduced_design,
          quiet = TRUE
        ) %>%
        results(
          tidy = TRUE
        ) %>%
        # Caution, results may contain NAs, seems to mainly depend on whether
        # zero cols were cleaned beforehand.
        return()
    },
    .is_uniform = function(x) {
      # Test whether all elements of x are the same.
      all(x == x[1])
    }
  )
)
