library("here")
library("R6")

suppressMessages(
  here::i_am("R6/Params.R")
)

Params <- R6Class(
  "Params"
)
