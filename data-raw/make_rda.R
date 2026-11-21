# data-raw/make_rda.R
# Produit data/ATcT.rda, data/RMG.rda, data/CHNOSZ.rda
# Les 3 sources restent SÉPARÉES (aucune fusion).

library(readr); library(dplyr)

dir.create("data", showWarnings = FALSE)

make_rda <- function(csv_path, prefix, obj_name, rda_path) {
  df <- read_csv(csv_path, show_col_types = FALSE) %>%
    mutate(
      compound_id = paste0(prefix, "_", row_number()),
      CAS         = as.character(CAS),
      delta_fH_0K = as.numeric(delta_fH_0K),
      uncertainty = as.numeric(uncertainty)
    )
  assign(obj_name, df)
  save(list = obj_name, file = rda_path, compress = "xz")
  message(rda_path, " : ", nrow(df), " lignes.")
}

make_rda(
  csv_path = "inst/extdata/ATcT/ATcT_consolidated.csv",
  prefix   = "ATcT",
  obj_name = "atct_harmonized",
  rda_path = "data/ATcT.rda"
)

make_rda(
  csv_path = "inst/extdata/RMG/RMG_consolidated.csv",
  prefix   = "RMG",
  obj_name = "rmg_harmonized",
  rda_path = "data/RMG.rda"
)

make_rda(
  csv_path = "inst/extdata/OBIGT/CHNOSZ_consolidated.csv",
  prefix   = "CHNOSZ",
  obj_name = "chnosz_harmonized",
  rda_path = "data/CHNOSZ.rda"
)