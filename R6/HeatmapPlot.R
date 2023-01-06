library("dplyr")
library("ggplot2")
library("here")
library("R6")

suppressMessages(
  here::i_am("R6/HeatmapPlot.R")
)

source(here("functions/norm_functions.R"))

source(here("R6/Plot.R"))

HeatmapPlot <- R6Class(
  "HeatmapPlot",
  inherit = Plot,
  public = list(
    initialize = function(reference_matrix,
                          title = NULL,
                          feature_scale = TRUE,
                          transcripts_ordered = NULL,
                          celltypes_ordered = NULL) {
      super$initialize(
        facet_width = 0.5
      )

      private$.check_matrix(reference_matrix)

      # TODO Add some checks on the inputs.
      private$.reference_matrix <- reference_matrix
      private$.title <- title
      private$.feature_scale <- feature_scale
      private$.transcript_order <- transcripts_ordered
      private$.celltype_order <- celltypes_ordered
    },
    save = function(dir, filename = "heatmap.png") {
      super$save(dir, filename)
    }
  ),
  active = list(
    reference_matrix = function() {
      private$.reference_matrix
    },
    title = function() {
      private$.title
    },
    feature_scale = function() {
      private$.feature_scale
    },
    plot_width = function() {
      self$params$base_width +
        self$params$facet_width * ncol(self$reference_matrix)
    },
    plot_height = function() {
      self$params$base_height + 20
    },
    count_df = function() {
      if (is.null(private$.count_df)) {
        private$.compute_count_df()
      }
      private$.count_df
    },
    plot = function() {
      if (is.null(private$.plot)) {
        private$.compute_plot()
      }
      private$.plot
    },
    transcript_order = function() {
      if (is.null(private$.transcript_order)) {
        private$.compute_transcript_order()
      }
      private$.transcript_order
    },
    celltype_order = function() {
      if (is.null(private$.celltype_order)) {
        private$.compute_celltype_order()
      }
      private$.celltype_order
    }
  ),
  private = list(
    .reference_matrix = NULL,
    .title = NULL,
    .feature_scale = NULL,
    .transcript_order = NULL,
    .celltype_order = NULL,
    .count_df = NULL,
    .plot = NULL,
    .compute_count_df = function() {
      count_df <- self$reference_matrix %>%
        t() %>%
        as.data.frame() %>%
        mutate(
          celltype = rownames(.)
        ) %>%
        select(celltype, everything()) %>%
        pivot_longer(
          !all_of(c("celltype")),
          names_to = "transcript",
          values_to = "abundance"
        )

      if (private$.feature_scale) {
        count_df <- count_df %>%
          group_by(transcript) %>%
          mutate(abundance = feature_scale(abundance))
      }

      private$.count_df <- count_df
    },
    .compute_celltype_order = function() {
      private$.celltype_order <- self$count_df %>%
          group_by(celltype) %>%
          summarize(
            mean = mean(abundance),
            .groups = "drop_last"
          ) %>%
          arrange(desc(mean)) %>%
          magrittr::extract2("celltype")
    },
    .compute_transcript_order = function() {
      # Transcripts will be ordered such that the first celltypes have the
      # transcripts with the highest abundances closest to the top of the plot.
      private$.transcript_order <- self$reference_matrix %>%
          # FIXME Do the sorting in a more straightforward way.
          as.data.frame() %>%
          arrange(across(all_of(self$celltype_order))) %>%
          rownames()
    },
    .compute_plot = function() {
      private$.plot <- self$count_df %>%
        mutate(
          celltype = factor(celltype, levels = self$celltype_order),
          transcript = factor(transcript, levels = self$transcript_order),
        ) %>%
        ggplot(
          aes(
            x = celltype,
            y = transcript,
            fill = abundance
          )
        ) +
        geom_raster() +
        labs(
          x = "Celltype",
          y = "Transcript",
          title = self$title
        ) +
        # TODO Make gradient adjustable.
        scale_fill_gradient(
          name = "Abundance",
          low = "dodgerblue4",
          high = "yellow"
        ) +
        self$themes$theme_benchmark() +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          ),
          axis.text.y = element_blank(),
          legend.text = element_text(
            angle = 45,
            hjust = 1
          ),
          legend.position = "bottom"
        )
    },
    .check_matrix = function(ref_matrix) {
      # ref_matrix cols are assumed to be celltypes, rows are assumed to be
      # transcripts
      stopifnot(
        !is.null(rownames(ref_matrix)),
        !is.null(colnames(ref_matrix)),
        !any(dim(ref_matrix) == 0)
      )
    }
  )
)