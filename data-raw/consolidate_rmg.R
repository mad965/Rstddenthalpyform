# data-raw/consolidate_rmg.R
source("data-raw/utils_harmonization.R")
library(readr); library(dplyr); library(stringr)

rmg_raw <- read_csv("../RMG/rmg_entries.csv", show_col_types = FALSE)

rmg_consolidated <- rmg_raw %>%
  transmute(
    compound_id        = NA_character_,
    name               = name,
    formula            = formula,
    CAS                = NA_character_,
    phase              = case_when(
      is.na(phase)        ~ NA_character_,
      phase == "gas"      ~ "g",
      phase == "liquid"   ~ "l",
      phase == "solid"    ~ "cr",
      phase == "aqueous"  ~ "aq",
      TRUE                ~ phase
    ),
    compound_nature    = ifelse(
      is.na(nature) | nature == "",
      classify_nature(formula, name),
      nature
    ),
    compound_subnature = classify_subnature(formula, name),
    delta_fH           = dfH_J_mol / 1000,
    delta_fG           = NA_real_,                       # non fourni
    S                  = as.numeric(S_J_mol_K),          # J/(mol·K), tel quel
    Cp                 = NA_real_,                       # non fourni
    delta_fH_0K        = NA_real_,
    temperature        = T_ref_K,
    uncertainty        = as.numeric(dfH_uncertainty_J_mol) / 1000,
    unit               = "kJ/mol",
    source             = "RMG",
    source_id          = id,
    reference          = reference
  ) %>%
  normalize_schema()

write_csv(rmg_consolidated, "../RMG/RMG_consolidated.csv")
message("RMG consolidé : ", nrow(rmg_consolidated), " lignes.")