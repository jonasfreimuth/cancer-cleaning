library("here")
library("ggplot2")

here::i_am("etc/plot_themes.R")

theme_benchmark <- function(...) {
  
  theme_minimal() %+replace%
    theme(panel.grid = element_blank(),
          axis.line = element_line(), ...)
  
}
