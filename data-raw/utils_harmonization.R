# data-raw/utils_harmonization.R
# Fonctions communes aux scripts de consolidation ATcT / RMG / CHNOSZ
# Aucune fusion entre sources : uniquement des helpers de nettoyage.

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
})

# --- Extraction de la phase depuis un texte type "CH4  (g)" ou "H2O  (cr,l)" ---
extract_phase <- function(x) {
  m <- str_match(x, "\\(([^)]+)\\)")
  phase <- m[, 2]
  dplyr::case_when(
    is.na(phase)             ~ NA_character_,
    str_detect(phase, "^g")  ~ "g",
    str_detect(phase, "^l")  ~ "l",
    str_detect(phase, "^cr") ~ "cr",
    str_detect(phase, "^aq") ~ "aq",
    TRUE                     ~ phase
  )
}

# --- Classification chimique (règles, pas "contient C -> organique") ---
classify_nature <- function(formula, name) {
  has_C     <- str_detect(formula, "C")
  has_H     <- str_detect(formula, "H")
  has_metal <- str_detect(formula,
                          "(Fe|Cu|Zn|Ni|Co|Mn|Cr|Pt|Pd|Ag|Au|Hg|Pb|Sn|Al|Mg|Ca|Na|K|Li)")
  
  dplyr::case_when(
    str_detect(formula, "CO3")                                           ~ "inorganic",
    formula %in% c("CO2", "CO")                                          ~ "inorganic",
    str_detect(name, "(?i)carbide")                                      ~ "inorganic",
    has_C & has_metal                                                    ~ "organometallic",
    str_detect(name, "(?i)(glycine|alanine|glucose|sucrose|urea|amino)") ~ "biochemical",
    has_C & has_H                                                        ~ "organic",
    TRUE                                                                 ~ "inorganic"
  )
}

classify_subnature <- function(formula, name) {
  dplyr::case_when(
    str_detect(name, "(?i)acid")                              ~ "acide carboxylique",
    str_detect(name, "(?i)alcohol|ol$")                       ~ "alcool",
    str_detect(name, "(?i)aldehyde")                          ~ "aldéhyde",
    str_detect(name, "(?i)ketone|one$")                       ~ "cétone",
    str_detect(name, "(?i)amine")                             ~ "amine",
    str_detect(name, "(?i)ester")                             ~ "ester",
    str_detect(name, "(?i)ether")                             ~ "éther",
    str_detect(name, "(?i)oxide")                             ~ "oxyde",
    str_detect(name, "(?i)sulfide")                           ~ "sulfure",
    str_detect(name, "(?i)chloride|bromide|iodide|fluoride")  ~ "halogénure",
    str_detect(name, "(?i)carbonate")                         ~ "carbonate",
    str_detect(name, "(?i)nitrate|sulfate|phosphate")         ~ "sel",
    TRUE                                                      ~ NA_character_
  )
}

# --- Schéma harmonisé commun (pour vérifier que les 3 sources ont les mêmes colonnes) ---
schema_harmonise <- c(
  "compound_id", "name", "formula", "CAS", "phase",
  "compound_nature", "compound_subnature",
  "delta_fH", "delta_fG", "S", "Cp",
  "delta_fH_0K", "temperature",
  "uncertainty", "unit", "source", "source_id", "reference"
)

# --- Ajoute les colonnes manquantes et réordonne ---
normalize_schema <- function(df) {
  for (col in schema_harmonise) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df[, schema_harmonise, drop = FALSE]
}

message("utils_harmonization.R chargé.")