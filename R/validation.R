#' Integrity checks of the database
#'
#' Runs structural checks on the normalised table: unique identifiers, valid
#' phases, reference conditions, formulas parsed, and consistency of the unit
#' conversion with the original values.
#'
#' @param db A database table, by default [thermo_db()].
#' @return A data frame with the columns `check` and `passed`.
#' @examples
#' check_database()
#' @export
check_database <- function(db = thermo_db()) {
  fac <- ifelse(db$orig_unit == "cal", .CAL_J, 1)
  checks <- list(
    ids_unique = !anyDuplicated(db$id),
    phases_valid = all(db$phase %in% c("aq", "cr", "liq", "gas")),
    reference_conditions = all(db$T_ref_K == 298.15 & db$p_ref_bar == 1),
    formulas_parsed = !anyNA(db$composition),
    dfH_conversion = isTRUE(all.equal(db$dfH_J_mol, db$orig_dfH * fac)),
    dfG_conversion = isTRUE(all.equal(db$dfG_J_mol, db$orig_dfG * fac)),
    has_dfH_values = any(!is.na(db$dfH_J_mol))
  )
  data.frame(check = names(checks), passed = unlist(checks),
             row.names = NULL, stringsAsFactors = FALSE)
}

# Entropy of the elements (J/(K mol of atoms)), including the pseudo-element
# Z of the electric charge, from element.csv of CHNOSZ. The CWM89 entries
# (Cox et al., 1989) are preferred.
.element_entropies <- function() {
  if (!is.null(.state$elements)) return(.state$elements)
  path <- file.path(.obigt_dir(), "element.csv")
  if (!file.exists(path)) {
    stop("Table des \u00e9l\u00e9ments introuvable.", call. = FALSE)
  }
  e <- utils::read.csv(path, colClasses = "character",
                       na.strings = c("NA", ""), stringsAsFactors = FALSE)
  e$s <- suppressWarnings(as.numeric(e$s))
  e$n <- suppressWarnings(as.numeric(e$n))
  e <- e[!is.na(e$s) & !is.na(e$n), , drop = FALSE]
  e <- e[order(e$source != "CWM89"), , drop = FALSE]
  e <- e[!duplicated(e$element), , drop = FALSE]
  s_atom <- stats::setNames(e$s / e$n * .CAL_J, e$element)
  .state$elements <- s_atom
  s_atom
}

.formula_entropy <- function(formula, s_atom) {
  p <- parse_formula(formula)
  nm <- names(p$elements)
  if (is.null(nm)) nm <- character(0)
  if (length(nm) && !all(nm %in% names(s_atom))) return(NA_real_)
  sum(p$elements * s_atom[nm]) + p$charge * s_atom[["Z"]]
}

#' Check the consistency of G, H and S
#'
#' For each entry that has the three values, compares the reported Gibbs
#' energy of formation with the one recomputed from the enthalpy of formation
#' and the entropy, at 298.15 K:
#' `dfG = dfH - T * (S - sum of the entropies of the elements)`. The entropies
#' of the elements (and of the pseudo-element of charge) come from the element
#' table of CHNOSZ. A large deviation does not prove an error but flags an
#' entry to examine (different conventions, inconsistent sources...).
#'
#' @param db A database table, by default [thermo_db()].
#' @param tol Tolerance in kJ/mol above which an entry is flagged.
#' @return A data frame sorted by decreasing absolute deviation, with
#'   `dfG_reported_kJ`, `dfG_recomputed_kJ`, `deviation_kJ` and `flagged`.
#'   Entries containing an element without entropy are not evaluated.
#' @examples
#' res <- check_gibbs_consistency()
#' head(res)
#' @export
check_gibbs_consistency <- function(db = thermo_db(), tol = 5) {
  s_atom <- .element_entropies()
  idx <- which(!is.na(db$dfG_J_mol) & !is.na(db$dfH_J_mol) &
                 !is.na(db$S_J_mol_K))
  s_el <- vapply(idx, function(i) .formula_entropy(db$formula[i], s_atom),
                 numeric(1))
  recomputed <- db$dfH_J_mol[idx] - 298.15 * (db$S_J_mol_K[idx] - s_el)
  dev <- (db$dfG_J_mol[idx] - recomputed) / 1000
  res <- data.frame(
    id = db$id[idx],
    name = db$name[idx],
    formula = db$formula[idx],
    phase = db$phase[idx],
    dfG_reported_kJ = db$dfG_J_mol[idx] / 1000,
    dfG_recomputed_kJ = recomputed / 1000,
    deviation_kJ = dev,
    flagged = !is.na(dev) & abs(dev) > tol,
    stringsAsFactors = FALSE
  )
  res <- res[order(-abs(res$deviation_kJ)), , drop = FALSE]
  rownames(res) <- NULL
  res
}
