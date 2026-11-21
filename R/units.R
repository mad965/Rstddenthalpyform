# Thermochemical calorie: 1 cal = 4.184 J (exact, by definition). This is the
# factor used by CHNOSZ for the species whose E_units column is "cal".
.CAL_J <- 4.184

.energy_factors <- c(
  "J/mol" = 1,
  "kJ/mol" = 1000,
  "cal/mol" = .CAL_J,
  "kcal/mol" = 1000 * .CAL_J
)

#' Convert an energy between units
#'
#' Converts molar energies using the thermochemical calorie
#' (1 cal = 4.184 J).
#'
#' @param x Numeric vector.
#' @param from,to Unit names among `"J/mol"`, `"kJ/mol"`, `"cal/mol"` and
#'   `"kcal/mol"`.
#' @return A numeric vector.
#' @examples
#' convert_energy(1, "kcal/mol", "kJ/mol")
#' @export
convert_energy <- function(x, from = "J/mol", to = "kJ/mol") {
  check_energy_unit(from)
  check_energy_unit(to)
  x * .energy_factors[[from]] / .energy_factors[[to]]
}

check_energy_unit <- function(unit) {
  ok <- is.character(unit) && length(unit) == 1L && !is.na(unit) &&
    unit %in% names(.energy_factors)
  if (!ok) {
    stop("Unit\u00e9 non reconnue. Unit\u00e9s accept\u00e9es : ",
         paste(names(.energy_factors), collapse = ", "), ".", call. = FALSE)
  }
  invisible(unit)
}

# "kJ/mol" -> "kJ/(mol K)" (unit of entropy and heat capacity)
entropy_unit <- function(units) {
  sub("/mol", "/(mol K)", units, fixed = TRUE)
}
