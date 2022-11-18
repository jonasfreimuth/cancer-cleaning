suppressMessages(here::i_am("functions/rmse.R"))
source(here::here("functions/err_fun_common.R"))

rmse <- function(pred, act, no.act.zeros = FALSE, na.rm = TRUE) {
  processed_input <- err_fun_common(pred, act, no.act.zeros, na.rm)

  pred <- processed_input[["pred"]]
  act <- processed_input[["act"]]

  return(sqrt(mean((pred - act)^2, na.rm = na.rm)))
}
