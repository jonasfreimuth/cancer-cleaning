library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Plot.R")
)

Plot <- R6Class(
  "Plot",
  public = list(
    initialize = function(base_width = 3,
                          base_height = 2,
                          facet_width = 5.5,
                          facet_height = 2) {
      # TODO Refactor params into own param class after it is clear how much
      # extra params Plot subclasses need.
      param_list <- list(
        base_width = base_width,
        base_height = base_height,
        facet_width = facet_width,
        facet_height = facet_height
      )

      private$.check_params(param_list)

      private$.params <- param_list
    },
    save = function(dir, filename) {
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
    params = function() {
      private$.params
    },
    themes = function() {
      private$.themes
    },
    plot_width = function() {
      if (is.null(private$.plot_width)) {
        private$.compute_plot_width()
      }
      private$.plot_width
    },
    plot_height = function() {
      if (is.null(private$.plot_height)) {
        private$.compute_plot_height()
      }
      private$.plot_height
    }
  ),
  private = list(
    .plot = NULL,
    .plot_width = NULL,
    .plot_height = NULL,
    .params = NULL,
    .themes = list(
      theme_benchmark = function(...) {
        theme_minimal() %+replace%
          theme(
            panel.grid = element_blank(),
            axis.line = element_line(), ...
          )
      }
    ),
    .compute_plot_width = function() {
      stop(paste(
        "No class implementation found, and no superclass implementation",
        "exists for plot_width."
      ))
    },
    .compute_plot_height = function() {
      stop(paste(
        "No class implementation found, and no superclass implementation",
        "exists for plot_height."
      ))
    },
    .check_params = function(param_list) {
      with(
        param_list,
        stopifnot(
          is.numeric(base_width),
          is.numeric(base_height),
          is.numeric(facet_width),
          is.numeric(facet_height)
        )
      )
    }
  )
)
