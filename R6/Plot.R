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
    }
  ),
  active = list(
    params = function() {
      private$.params
    },
    themes = function() {
      private$.themes
    }
  ),
  private = list(
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
