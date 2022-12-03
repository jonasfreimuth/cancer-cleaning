library("data.table")
library("here")

here::i_am("scripts/generate_test_data.R")

source(here("functions/util_functions.R"))
source(here("functions/benchmark_functions.R"))

seed <- 123

set.seed(seed)

source_data_dir <- here("datasets/Wu_etal_2021_BRCA_scRNASeq/")
destin_data_dir <- here("datasets/Wu_etal_downsampled_test")

data <- load_experiment(
  count_mat_file = here(source_data_dir, "count_matrix_sparse.mtx"),
  rowname_file = here(source_data_dir, "count_matrix_genes.tsv"),
  colname_file = here(source_data_dir, "count_matrix_barcodes.tsv"),
  meta_file = here(source_data_dir, "metadata.csv"),
  testing = TRUE
)

dir.create(destin_data_dir, showWarnings = FALSE, recursive = TRUE)

writeMM(data$count_mat, here(destin_data_dir, "count_matrix_sparse.mtx"))

writeLines(
  rownames(data$count_mat),
  here(destin_data_dir, "count_matrix_genes.tsv")
)

writeLines(
  colnames(data$count_mat),
  here(destin_data_dir, "count_matrix_barcodes.tsv")
)

fwrite(data$meta, here(destin_data_dir, "metadata.csv"))
