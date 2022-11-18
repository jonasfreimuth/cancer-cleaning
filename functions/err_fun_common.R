suppressMessages(here::i_am("functions/err_fun_common.R"))

err_fun_common <- function(pred, act, no.act.zeros = FALSE, na.rm = TRUE) {
  if (no.act.zeros) {
    sel_vec <- act != 0
  } else {
    sel_vec <- TRUE
  }

  pred <- pred[sel_vec]
  act <- act[sel_vec]

  if (na.rm) {
    na_ind_vec <- is.na(pred) | is.na(act)
    pred <- pred[!na_ind_vec]
    act <- act[!na_ind_vec]
  }

  if (length(pred) != length(act)) {
    stop("pred and act must have the same length!")
  }

  return(list("pred" = pred, "act" = act))
}
