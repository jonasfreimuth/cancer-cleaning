suppressMessages(here::i_am("functions/kld.R"))
source(here::here("functions/err_fun_common.R"))

kld <- function(pred, act, no.act.zeros = FALSE, na.rm = TRUE) {
  processed_input <- err_fun_common(pred, act, no.act.zeros, na.rm)

  pred <- processed_input[["pred"]]
  act <- processed_input[["act"]]

  return(sum(pred * log(pred / act)))
}
