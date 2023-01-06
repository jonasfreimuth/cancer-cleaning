library("data.table")
library("dplyr")
library("ggplot2")
library("here")
library("R6")

suppressMessages(
  here::i_am("R6/RmsePlot.R")
)

source(here("R6/Plot.R"))

RmsePlot <- R6Class(
  "RmsePlot",
  inherit = Plot,
  public = list(
    initialize = function(celltype_rmse_df, summary_rmse_df) {
      # Unsure if this is necessary.
      # TODO Add global plotting params (Currently this is done via the default)
      # params to Plot
      super$initialize()

      private$.check_celltype_rmse_df(celltype_rmse_df)
      private$.check_summary_rmse_df(summary_rmse_df)

      private$.celltype_rmse_df <- celltype_rmse_df
      private$.summary_rmse_df <- summary_rmse_df
    },
    save = function(dir, filename = "rmse_plot.png") {
      path <- here(dir, filename)

      self$plot %>%
        {
          ggsave(
            plot = .,
            filename = path,
            width = self$plot_width,
            height = self$plot_height,
            limitsize = FALSE
          )
        }
    }
  ),
  active = list(
    plot_width = function() {
      col_df <- self$plot_df %>%
        select(all_of(self$facet_list[["cols"]]))

      n_plot_cols <- apply(
        col_df,
        2,
        data.table::uniqueN
      ) %>%
        sum()

      self$params$base_width + self$params$facet_width * n_plot_cols
    },
    plot_height = function() {
      row_df <- self$plot_df %>%
        select(all_of(self$facet_list$rows))

      n_plot_rows <- apply(
        row_df,
        2,
        data.table::uniqueN
      ) %>%
        sum()

      self$params$base_height + self$params$facet_height * n_plot_rows
    },
    celltype_rmse_df = function() {
      private$.celltype_rmse_df
    },
    summary_rmse_df = function() {
      private$.summary_rmse_df
    },
    plot_df = function() {
      if (is.null(private$.plot_df)) {
        private$.compute_plot_df()
      }
      private$.plot_df
    },
    facet_params = function() {
      if (is.null(private$.facet_params)) {
        private$.compute_facet_params()
      }
      private$.facet_params
    },
    facet_list = function() {
      if (is.null(private$.facet_list)) {
        private$.compute_facet_list()
      }
      private$.facet_list
    },
    facet_formula = function() {
      if (is.null(private$.facet_formula)) {
        private$.compute_facet_formula()
      }
      private$.facet_formula
    },
    plot = function() {
      if (is.null(private$.plot)) {
        private$.compute_plot()
      }
      private$.plot
    }
  ),
  private = list(
    .celltype_rmse_df = NULL,
    .summary_rmse_df = NULL,
    .plot_df = NULL,
    .facet_params = NULL,
    .facet_list = NULL,
    .facet_formula = NULL,
    .plot = NULL,
    .compute_plot_df = function() {
      private$.plot_df <- self$celltype_rmse_df %>%
        left_join(
          self$summary_rmse_df,
          by = self$facet_params
        ) %>%
        ungroup()
    },
    .compute_facet_params = function() {
      facet_params <- intersect(
        names(self$celltype_rmse_df),
        names(self$summary_rmse_df)
      )

      private$.check_facet_params(facet_params)

      private$.facet_params <- facet_params
    },
    .compute_facet_list = function() {
      n_facet_params <- length(self$facet_params)
      cutoff <- ceiling(n_facet_params / 2)

      private$.facet_list <- list(
        cols = self$facet_params[1:cutoff],
        rows = self$facet_params[(cutoff + 1):n_facet_params]
      )
    },
    .compute_facet_formula = function() {
      private$.facet_formula <- as.formula(paste0(
        paste(self$facet_list$cols, collapse = " + "),
        " ~ ",
        paste(self$facet_list$rows, collapse = " + ")
      ))
    },
    .compute_plot = function() {
      private$.plot <- self$plot_df %>%
        ggplot(
          aes(celltype, rmse)
        ) +
        geom_col(alpha = 0.5, position = position_identity()) +
        geom_text(
          aes(
            x = length(unique(celltype)) / 2,
            y = max(rmse) * 1.1,
            label = paste0(
              "Mean RMSE: ", round(mean_rmse, 3)
            )
          ),
          vjust = 1
        ) +
        labs(
          x = "Cell types",
          y = "Per celltype RMSE"
        ) +
        facet_grid(
          # FIXME Figure out a way to include labels for each var.
          rows = self$facet_formula
        ) +
        self$themes$theme_benchmark() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    },
    .check_celltype_rmse_df = function(celltype_rmse_df) {
      stopifnot(
        is.data.frame(celltype_rmse_df),
        c("celltype", "rmse") %in% names(celltype_rmse_df),
        nrow(celltype_rmse_df) >= 1,
        is.numeric(celltype_rmse_df$rmse),
        all(
          sapply(select(celltype_rmse_df, -rmse), class) %in%
            c("character", "factor")
        )
      )
    },
    .check_summary_rmse_df = function(summary_rmse_df) {
      stopifnot(
        is.data.frame(summary_rmse_df),
        c("mean_rmse") %in% names(summary_rmse_df),
        nrow(summary_rmse_df) >= 1,
        is.numeric(summary_rmse_df$mean_rmse),
        all(
          sapply(select(summary_rmse_df, -mean_rmse), class) %in%
            c("character", "factor")
        )
      )
    },
    .check_facet_params = function(facet_params) {
      stopifnot(
        length(facet_params) >= 1
      )
    }
  )
)
