suppressMessages(here::i_am("functions/med_abs_err.R"))
source(here::here("functions/err_fun_common.R"))

med_abs_err <- function(pred, act, no.act.zeros = FALSE, na.rm = TRUE) {
  processed_input <- err_fun_common(pred, act, no.act.zeros, na.rm)

  pred <- processed_input[["pred"]]
  act <- processed_input[["act"]]

  return(median(abs(pred - act), na.rm = na.rm))
}
