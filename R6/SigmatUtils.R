library("DESeq2")
library("scuttle")
library("here")
library("R6")

here::i_am("R6/SigmatUtils.R")

SigmatUtils <- R6Class(
  "SigmatUtils",
  public = list(
    create_marker_df = function(matrix, meta) {
      matrix %>%
        # remove rows that are all the same
        # Also includes rows that are all zero-counts
        magrittr::extract(
          i = !apply(
            ., 1, private$is_uniform
          ), j = , drop = FALSE
        ) %>%
        private$get_de_transcripts(
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
    }
  ),
  private = list(
    get_de_transcripts = function(count_mat, meta, design) {
      # TODO Get the functions in here to STFU.
      meta <- meta %>%
        arrange(cell) %>%
        # Prevent DESeq fun from complaining
        mutate(across(where(is.character), as.factor))
      count_mat <- count_mat[, order(colnames(count_mat))]

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
    is_uniform = function(x) {
      # Test whether all elements of x are the same.
      all(x == x[1])
    }
  )
)
