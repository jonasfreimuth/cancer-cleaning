library("here")
library("R6")

here::i_am("R6/ScRnaExperiment.R")

source(here("R6/CountMatrix"))

ScRnaExperiment <- R6Class(
  "ScRnaExperiment",
  inherit = CountMatrix,
  public(
    initialize = function(count_mat_file, rowname_file, colname_file,
                          meta_file,
                          cell_col = "V1", celltype_col = "celltype_major",
                          downsample_frac) {
      rename_vec <- c(
        "cell" = cell_col,
        "celltype" = celltype_col
      )

      self$meta <- fread(meta_file) %>%
        # See rlang::`!!!`.
        rename(!!!rename_vec) %>%
        # Ordering ensured here.
        arrange(cell)

      self$matrix <- readMM(count_mat_file)

      transcripts <- readLines(rowname_file)
      cells <- readLines(colname_file)

      self$matrix@Dimnames <- list(
        transcripts,
        cells
      )

      # Ordering ensured here.
      self$matrix <- self$matrix[order(self$cells), ]

      if (!base::missing(downsample_frac)) {
        private$downsample()
      }
    }
  ),
  private = list(
    downsample = function(cell_frac, transcript_frac) {
      rnd_cell_idx <- self$n_cells %>%
        {
          sample(x = seq_len(.), size = cell_frac * .)
        }
      rnd_transcript_idx <- self$n_celltypes %>%
        {
          sample(x = seq_len(.), size = transcript_frac * .)
        }

      self$matrix <- self$matrix[rnd_cell_idx, rnd_transcript_idx]

      self$meta <- self$meta %>%
        filter(cell %in% self$cells)
    }
  )
)
