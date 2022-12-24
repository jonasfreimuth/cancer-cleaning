library("styler")

list_dir <- function(path) {
  dir(path, full.names = TRUE, pattern = "*.R")
}

style_file(c(
  list_dir("functions"),
  list_dir("scripts"),
  list_dir("aux_scripts"),
  list_dir("R6")
))
