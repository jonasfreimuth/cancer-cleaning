library("dplyr")
library("ggplot2")
library("here")
library("R6")

suppressMessages(
  here::i_am("R6/PlotUtils.R")
)

source(here("functions/norm_functions.R"))
source(here("etc/plot_themes.R"))

PlotUtils <- R6Class(
  "PlotUtils",
  public = list(
    create_heatmap = function(reference_matrix, title = NULL,
                              feature_scale = TRUE,
                              transcripts_ordered = NULL,
                              celltypes_ordered = NULL) {
      count_df <- private$.count_df_from_ref_matrix(
        reference_matrix,
        feature_scale
      )

      if (is.null(celltypes_ordered)) {
        celltypes_ordered <- count_df %>%
          group_by(celltype) %>%
          summarize(
            mean = mean(abundance),
            .groups = "drop_last"
          ) %>%
          arrange(desc(mean)) %>%
          magrittr::extract2("celltype")
      }


      if (is.null(transcripts_ordered)) {
        transcripts_ordered <- reference_matrix %>%
          # FIXME Do the sorting in a  more straightforward way.
          as.data.frame() %>%
          arrange(across(all_of(celltypes_ordered))) %>%
          rownames()
      }

      qc_plot <- count_df %>%
        mutate(
          celltype = factor(celltype, levels = celltypes_ordered),
          transcript = factor(transcript, levels = transcripts_ordered),
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
          title = title
        ) +
        scale_fill_gradient(
          name = "Abundance",
          low = "dodgerblue4",
          high = "yellow"
        ) +
        theme_benchmark() +
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

      return(qc_plot)
    }
  ),
  private = list(
    .count_df_from_ref_matrix = function(ref_matrix, feature_scale = FALSE) {
      # ref_matrix cols are assumed to be celltypes, rows are assumed to be
      # transcripts
      stopifnot(
        !is.null(rownames(ref_matrix)),
        !is.null(colnames(ref_matrix)),
        !any(dim(ref_matrix) == 0)
      )

      count_df <- ref_matrix %>%
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

      if (feature_scale) {
        count_df <- count_df %>%
          group_by(transcript) %>%
          mutate(abundance = feature_scale(abundance))
      }

      return(count_df)
    }
  )
)
