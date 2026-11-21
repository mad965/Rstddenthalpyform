#' Create a predictor
#'
#' A predictor holds the database used by the lookup functions. Creating one
#' makes it the default predictor used by the other functions.
#'
#' @param db Database table, by default [thermo_db()].
#' @return An object of class `ThermoPredictor`.
#' @examples
#' thermo <- ThermoPredictor()
#' thermo
#' @export
ThermoPredictor <- function(db = thermo_db()) {
  obj <- structure(list(db = db), class = "ThermoPredictor")
  .state$predictor <- obj
  obj
}

#' @export
print.ThermoPredictor <- function(x, ...) {
  cat("<ThermoPredictor>\n")
  cat("  Base de donn\u00e9es : ", nrow(x$db), " entr\u00e9es (",
      unique(x$db$source)[1L], ")\n", sep = "")
  invisible(x)
}

default_predictor <- function() {
  if (is.null(.state$predictor)) ThermoPredictor()
  .state$predictor
}

.bind_results <- function(lst) {
  res <- do.call(rbind, lapply(lst, function(r) {
    class(r) <- "data.frame"
    r
  }))
  rownames(res) <- NULL
  class(res) <- c("thermo_result", "data.frame")
  res
}

.formation_database <- function(sp, phase, units, db) {
  sel <- .select_rows(sp, phase, db)
  idx <- sel$index
  if (!length(idx)) {
    q <- trimws(sp)
    hits <- tryCatch(agrep(q, db$name, max.distance = 0.15, ignore.case = TRUE),
                     error = function(e) integer(0))
    sugg <- unique(db$name[hits])
    note <- "Esp\u00e8ce absente de la base"
    if (!is.null(phase)) note <- paste0(note, " pour la phase demand\u00e9e")
    note <- paste0(note, ".")
    if (length(sugg)) {
      note <- paste0(note, " Suggestions : ",
                     paste(utils::head(sugg, 5L), collapse = ", "), ".")
    }
    return(new_thermo_result(query = sp, unit = units, method = "database",
                             method_class = "database",
                             source = unique(db$source)[1L],
                             status = "not_found", note = note))
  }
  fac <- .energy_factors[[units]]
  v <- db$dfH_J_mol[idx] / fac
  ref <- ifelse(is.na(db$ref2[idx]), db$ref1[idx],
                paste(db$ref1[idx], db$ref2[idx], sep = "; "))
  note <- rep("", length(idx))
  note[is.na(v)] <- paste0(
    "Esp\u00e8ce pr\u00e9sente mais sans valeur de \u0394fH\u00b0 dans la base ",
    "(la valeur n'est pas invent\u00e9e : elle reste manquante).")
  extra <- character(0)
  if (identical(sel$match, "composition")) {
    extra <- c(extra, paste0("Correspondance par composition \u00e9l\u00e9mentaire : ",
                             "la formule saisie diff\u00e8re de celle de la base."))
  }
  n_names <- tapply(db$name[idx], db$phase[idx],
                    function(z) length(unique(z)))
  if (any(n_names > 1L)) {
    extra <- c(extra, paste0("Plusieurs esp\u00e8ces (isom\u00e8res ?) correspondent \u00e0 ",
                             "la requ\u00eate pour une m\u00eame phase : v\u00e9rifier la ",
                             "colonne name."))
  }
  if (length(extra)) {
    note <- ifelse(nzchar(note), paste(note, paste(extra, collapse = " ")),
                   paste(extra, collapse = " "))
  }
  new_thermo_result(
    query = sp, id = db$id[idx], name = db$name[idx],
    formula = db$formula[idx], nature = db$nature[idx], phase = db$phase[idx],
    value = v, unit = units, T_K = db$T_ref_K[idx], p_bar = db$p_ref_bar[idx],
    method = "database", method_class = "database", source = db$source[idx],
    reference = ref, uncertainty = NA_real_,
    uncertainty_note = "Non fournie par OBIGT",
    status = ifelse(is.na(v), "no_value", "ok"), note = note
  )
}

#' Standard enthalpy of formation of a species
#'
#' Looks for the species in the database and returns the documented value of
#' its standard enthalpy of formation with its source.
#'
#' The result always gives the value, the unit, the phase, the reference
#' temperature, the method, the source and the uncertainty when available.
#' OBIGT does not provide uncertainties: the column `uncertainty` is `NA`.
#' When several entries match (several phases, isomers), all are returned: the
#' phase and the reference conditions are never guessed. A species that is not
#' in the database gets the status `"not_found"`, and a species present
#' without a value gets the status `"no_value"`: no value is ever invented.
#'
#' @param species Character vector of names or formulas.
#' @param phase Optional phase (`"aq"`, `"cr"`, `"liq"`, `"gas"`).
#' @param units Energy unit: `"kJ/mol"` (default), `"J/mol"`, `"kcal/mol"` or
#'   `"cal/mol"`.
#' @param predictor A [ThermoPredictor()]; by default the current one.
#' @return A data frame of class `thermo_result`, see [new_thermo_result()].
#' @examples
#' get_formation_enthalpy("methane", phase = "gas")
#' get_formation_enthalpy("CH4", phase = "gas", units = "kcal/mol")
#' @export
get_formation_enthalpy <- function(species, phase = NULL, units = "kJ/mol",
                                   predictor = default_predictor()) {
  check_energy_unit(units)
  phase <- .normalize_phase(phase)
  if (!is.character(species) || !length(species)) {
    stop("`species` doit \u00eatre un vecteur de cha\u00eenes de caract\u00e8res.",
         call. = FALSE)
  }
  .bind_results(lapply(species, .formation_database, phase = phase,
                       units = units, db = predictor$db))
}
