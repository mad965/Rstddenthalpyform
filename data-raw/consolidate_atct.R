# data-raw/consolidate_atct.R
source("data-raw/utils_harmonization.R")
library(readr); library(dplyr); library(stringr)

atct_raw <- read_csv("../ATcT/ATcT_enthalpies.csv", show_col_types = FALSE)

atct_consolidated <- atct_raw %>%
  mutate(
    phase         = extract_phase(formule),
    formula_clean = str_trim(str_remove(formule, "\\s*\\([^)]*\\)")),
    name_clean    = str_trim(str_remove(nom,    "\\s*\\([^)]*\\)")),
    CAS_clean     = str_remove(CAS, "^\\*")
  ) %>%
  transmute(
    compound_id        = NA_character_,
    name               = name_clean,
    formula            = formula_clean,
    CAS                = CAS_clean,
    phase              = phase,
    compound_nature    = classify_nature(formula_clean, name_clean),
    compound_subnature = classify_subnature(formula_clean, name_clean),
    delta_fH           = Hf_298K,
    delta_fG           = NA_real_,   # non fourni par ATcT
    S                  = NA_real_,   # non fourni par ATcT
    Cp                 = NA_real_,   # non fourni par ATcT
    delta_fH_0K        = Hf_0K,
    temperature        = 298.15,
    uncertainty        = incertitude_298K,
    unit               = "kJ/mol",
    source             = "ATcT",
    source_id          = id,
    reference          = paste0("ATcT TN ", version)
  ) %>%
  normalize_schema()

write_csv(atct_consolidated, "../ATcT/ATcT_consolidated.csv")
message("ATcT consolidé : ", nrow(atct_consolidated), " lignes.")