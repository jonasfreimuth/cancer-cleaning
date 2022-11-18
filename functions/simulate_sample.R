library("magrittr")
library("dplyr")


simulate_sample <- function(n_vars, n_muts) {
  vars <- LETTERS[1:n_vars]
  muts <- letters[1:n_muts]

  true_vals <- data.frame(
    variant = vars,
    prop = runif(n_vars)
  ) %>%
    mutate(prop = prop / sum(prop))

  sigmut_mat <- matrix(
    sample(c(0, 1), n_muts * n_vars, replace = TRUE),
    nrow = n_muts,
    ncol = n_vars,
    dimnames = list(
      muts, vars
    )
  )

  mut_prop <- sigmut_mat %*% as.matrix(true_vals$prop)

  return(list(sigmut_mat, mut_prop, true_vals))
}
