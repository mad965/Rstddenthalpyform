# data-raw/consolidate_chnosz.R
source("data-raw/utils_harmonization.R")
library(readr); library(dplyr); library(CHNOSZ)

obigt <- import_obigt(dir = "inst/extdata/OBIGT")

E.units("J")

obigt$G_subcrt  <- NA_real_
obigt$H_subcrt  <- NA_real_
obigt$S_subcrt  <- NA_real_
obigt$Cp_subcrt <- NA_real_

for (i in seq_len(nrow(obigt))) {
  sp <- obigt$name[i]
  st <- obigt$state[i]
  out <- tryCatch(
    subcrt(sp, state = st, T = 25, P = 1,
           property = c("G", "H", "S", "Cp")),
    error = function(e) NULL
  )
  if (!is.null(out) && length(out$out) > 0) {
    o <- out$out[[1]]
    obigt$G_subcrt[i]  <- if (!is.null(o$G))  o$G[1]  else NA_real_
    obigt$H_subcrt[i]  <- if (!is.null(o$H))  o$H[1]  else NA_real_
    obigt$S_subcrt[i]  <- if (!is.null(o$S))  o$S[1]  else NA_real_
    obigt$Cp_subcrt[i] <- if (!is.null(o$Cp)) o$Cp[1] else NA_real_
  }
  if (i %% 500 == 0) message("  ", i, " / ", nrow(obigt))
}

chnosz_consolidated <- obigt %>%
  transmute(
    compound_id        = NA_character_,
    name               = name,
    formula            = formula,
    CAS                = NA_character_,
    phase              = case_when(
      phase == "gas" ~ "g",
      phase == "liq" ~ "l",
      TRUE ~ phase
    ),
    compound_nature    = nature,
    compound_subnature = classify_subnature(formula, name),
    delta_fH           = H_subcrt / 1000,
    delta_fG           = G_subcrt / 1000,
    S                  = S_subcrt,
    Cp                 = Cp_subcrt,
    delta_fH_0K        = NA_real_,
    temperature        = 298.15,
    uncertainty        = NA_real_,
    unit               = "kJ/mol",
    source             = "CHNOSZ",
    source_id          = id,
    reference          = ifelse(is.na(ref1), NA_character_, ref1)
  ) %>%
  normalize_schema()

write_csv(chnosz_consolidated, "inst/extdata/OBIGT/CHNOSZ_consolidated.csv")
message("CHNOSZ consolidé : ", nrow(chnosz_consolidated), " lignes, ",
        sum(is.na(chnosz_consolidated$delta_fH)), " sans delta_fH.")