#' Build a formation-enthalpy result table
#'
#' Constructor of the result object returned by [get_formation_enthalpy()].
#' It is exported so that future estimation methods can return results in the
#' same format, with `method_class = "estimate"`.
#'
#' All arguments are recycled to the length of the longest one.
#'
#' @param query Query typed by the user.
#' @param id,name,formula,nature,phase Identification of the species.
#' @param value Standard enthalpy of formation, in `unit`.
#' @param unit Unit of `value`.
#' @param T_K,p_bar Reference temperature (K) and pressure (bar).
#' @param method Name of the method (`"database"`, `"group_contribution"`...).
#' @param method_class `"database"` for a documented value, `"estimate"` for
#'   an estimate. An estimate is never an experimental value.
#' @param source,reference Source of the value and its bibliographic keys.
#' @param uncertainty,uncertainty_note Uncertainty (same unit as `value`) and a
#'   comment on it.
#' @param status `"ok"` (documented value), `"no_value"` (species present
#'   without a value) or `"not_found"` (species absent from the database).
#' @param note Free-text comment.
#' @return A data frame of class `thermo_result`.
#' @export
new_thermo_result <- function(query, id = NA_character_, name = NA_character_,
                              formula = NA_character_, nature = NA_character_,
                              phase = NA_character_, value = NA_real_,
                              unit = "kJ/mol", T_K = NA_real_, p_bar = NA_real_,
                              method = NA_character_,
                              method_class = NA_character_,
                              source = NA_character_, reference = NA_character_,
                              uncertainty = NA_real_,
                              uncertainty_note = NA_character_,
                              status = "ok", note = "") {
  n <- max(length(query), length(id), length(name), length(value),
           length(status), length(phase), length(method))
  df <- data.frame(
    query = rep_len(as.character(query), n),
    id = rep_len(as.character(id), n),
    name = rep_len(as.character(name), n),
    formula = rep_len(as.character(formula), n),
    nature = rep_len(as.character(nature), n),
    phase = rep_len(as.character(phase), n),
    value = rep_len(as.numeric(value), n),
    unit = rep_len(as.character(unit), n),
    T_K = rep_len(as.numeric(T_K), n),
    p_bar = rep_len(as.numeric(p_bar), n),
    method = rep_len(as.character(method), n),
    method_class = rep_len(as.character(method_class), n),
    source = rep_len(as.character(source), n),
    reference = rep_len(as.character(reference), n),
    uncertainty = rep_len(as.numeric(uncertainty), n),
    uncertainty_note = rep_len(as.character(uncertainty_note), n),
    status = rep_len(as.character(status), n),
    note = rep_len(as.character(note), n),
    stringsAsFactors = FALSE
  )
  class(df) <- c("thermo_result", "data.frame")
  df
}

#' @export
print.thermo_result <- function(x, ...) {
  df <- x
  class(df) <- "data.frame"
  if (nrow(df) == 0L) {
    cat("<thermo_result> : aucun r\u00e9sultat\n")
    return(invisible(x))
  }
  dash <- function(v) ifelse(is.na(v), "-", v)
  show <- data.frame(
    query = dash(df$query),
    name = dash(df$name),
    formula = dash(df$formula),
    phase = dash(df$phase),
    value = formatC(df$value, format = "f", digits = 2),
    unit = dash(df$unit),
    T_K = formatC(df$T_K, format = "f", digits = 2),
    method = dash(df$method),
    source = dash(df$source),
    uncertainty = ifelse(is.na(df$uncertainty), "n/a",
                         formatC(df$uncertainty, format = "f", digits = 2)),
    status = dash(df$status),
    stringsAsFactors = FALSE
  )
  print(show, row.names = FALSE, right = FALSE)
  notes <- unique(df$note[!is.na(df$note) & nzchar(df$note)])
  for (nt in notes) cat("Note : ", nt, "\n", sep = "")
  if (any(df$method_class %in% "database" & df$status %in% "ok")) {
    cat("Valeur document\u00e9e de la base de donn\u00e9es (non une pr\u00e9diction).",
        "Incertitude : n/a = non fournie par la source.\n")
  }
  if (any(df$method_class %in% "estimate" & df$status %in% "ok")) {
    cat("ESTIMATION : \u00e0 ne jamais pr\u00e9senter comme une valeur exp\u00e9rimentale.\n")
  }
  invisible(x)
}

# Rows of the database matching a query, optionally restricted to phases
.select_rows <- function(query, phase, db) {
  m <- find_species(query, db)
  idx <- m$index
  if (length(idx) && !is.null(phase)) idx <- idx[db$phase[idx] %in% phase]
  list(index = idx, match = m$match)
}

#' Standard thermodynamic properties of a species
#'
#' Returns, for each matching entry of the database, the standard enthalpy
#' and Gibbs energy of formation, the standard entropy, the heat capacity and
#' the volume at 298.15 K and 1 bar, with their source. The entropy and the
#' heat capacity are absolute standard values, not values of formation.
#' Missing values are `NA`; they are never replaced by zero.
#'
#' @param species Character vector of names or formulas.
#' @param phase Optional phase (`"aq"`, `"cr"`, `"liq"`, `"gas"`).
#' @param units Energy unit: `"kJ/mol"` (default), `"J/mol"`, `"kcal/mol"` or
#'   `"cal/mol"`. Entropy and heat capacity use the same energy unit per kelvin.
#' @param predictor A [ThermoPredictor()]; by default the current one.
#' @return A data frame (`NULL` with a warning if nothing is found) with the
#'   columns `dfH`, `dfG`, `S`, `Cp`, `V_cm3_mol` and their units.
#' @examples
#' get_thermo_properties("methane", phase = "gas")
#' @export
get_thermo_properties <- function(species, phase = NULL, units = "kJ/mol",
                                  predictor = default_predictor()) {
  check_energy_unit(units)
  phase <- .normalize_phase(phase)
  db <- predictor$db
  f <- .energy_factors[[units]]
  out <- list()
  for (sp in species) {
    sel <- .select_rows(sp, phase, db)
    i <- sel$index
    if (!length(i)) {
      warning("Aucune esp\u00e8ce trouv\u00e9e pour '", sp, "'.", call. = FALSE)
      next
    }
    out[[length(out) + 1L]] <- data.frame(
      query = sp,
      match = sel$match,
      id = db$id[i],
      name = db$name[i],
      formula = db$formula[i],
      nature = db$nature[i],
      phase = db$phase[i],
      state = db$state[i],
      dfH = db$dfH_J_mol[i] / f,
      dfG = db$dfG_J_mol[i] / f,
      S = db$S_J_mol_K[i] / f,
      Cp = db$Cp_J_mol_K[i] / f,
      V_cm3_mol = db$V_cm3_mol[i],
      energy_unit = units,
      entropy_unit = entropy_unit(units),
      T_K = db$T_ref_K[i],
      p_bar = db$p_ref_bar[i],
      source = db$source[i],
      ref1 = db$ref1[i],
      ref2 = db$ref2[i],
      stringsAsFactors = FALSE
    )
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
