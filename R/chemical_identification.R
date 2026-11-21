#' Parse a chemical formula
#'
#' Splits a formula written in the OBIGT/CHNOSZ notation into element counts
#' and electric charge. Supported syntax: element symbols with integer or
#' decimal subscripts, parentheses with multipliers, hydrates written with
#' `*` or `.nH2O`, adducts written with `:`, and a trailing charge such as
#' `+`, `-`, `+2` or `-3`. The formula `(Z-1)` of the electron is recognised.
#'
#' @param formula A single character string.
#' @return A list with `elements` (named numeric vector of atom counts) and
#'   `charge` (numeric).
#' @examples
#' parse_formula("CO3-2")
#' parse_formula("Mg2(OH)2(CO3)*3H2O")
#' @export
parse_formula <- function(formula) {
  if (!is.character(formula) || length(formula) != 1L || is.na(formula)) {
    stop("`formula` doit \u00eatre une cha\u00eene de caract\u00e8res unique.", call. = FALSE)
  }
  if (formula == "(Z-1)") {
    return(list(elements = stats::setNames(numeric(0), character(0)),
                charge = -1))
  }
  f <- formula
  charge <- 0
  m <- regexpr("[+-][0-9]*$", f)
  if (m > 0L) {
    ch <- substr(f, m, nchar(f))
    sign <- if (substr(ch, 1L, 1L) == "+") 1 else -1
    digits <- substr(ch, 2L, nchar(ch))
    charge <- sign * (if (nzchar(digits)) as.numeric(digits) else 1)
    f <- substr(f, 1L, m - 1L)
  }
  total <- stats::setNames(numeric(0), character(0))
  # hydrate written with a dot, e.g. "MgSO4.H2O" or "Yb2(CO3)3.6H2O"
  mh <- regexpr("\\.[0-9]*H2O$", f)
  if (mh > 0L) {
    ntxt <- sub("^\\.([0-9]*)H2O$", "\\1", substr(f, mh, nchar(f)))
    n <- if (nzchar(ntxt)) as.numeric(ntxt) else 1
    total <- .add_counts(total, c(H = 2 * n, O = n))
    f <- substr(f, 1L, mh - 1L)
  }
  # hydrates ("*") and adducts (":") are summed
  parts <- strsplit(f, "[*:]")[[1L]]
  for (part in parts) {
    coef <- 1
    mc <- regexpr("^[0-9]+", part)
    if (mc > 0L) {
      len <- attr(mc, "match.length")
      coef <- as.numeric(substr(part, 1L, len))
      part <- substr(part, len + 1L, nchar(part))
    }
    total <- .add_counts(total, .parse_group(part) * coef)
  }
  list(elements = total, charge = charge)
}

.tokenize_formula <- function(s) {
  m <- gregexpr("[A-Z][a-z]?|[0-9]+(\\.[0-9]+)?|\\(|\\)", s)
  tokens <- regmatches(s, m)[[1L]]
  if (!identical(paste(tokens, collapse = ""), s)) {
    stop("Formule non analysable : ", s, call. = FALSE)
  }
  tokens
}

.add_counts <- function(total, add) {
  for (el in names(add)) {
    if (el %in% names(total)) {
      total[[el]] <- total[[el]] + add[[el]]
    } else {
      total[el] <- add[[el]]
    }
  }
  total
}

.parse_group <- function(s) {
  tokens <- .tokenize_formula(s)
  n <- length(tokens)
  stack <- list(numeric(0))
  i <- 1L
  while (i <= n) {
    tk <- tokens[i]
    if (tk == "(") {
      stack[[length(stack) + 1L]] <- numeric(0)
      i <- i + 1L
    } else if (tk == ")") {
      i <- i + 1L
      mult <- 1
      if (i <= n && grepl("^[0-9]", tokens[i])) {
        mult <- as.numeric(tokens[i])
        i <- i + 1L
      }
      if (length(stack) < 2L) {
        stop("Parenth\u00e8ses non appari\u00e9es dans la formule : ", s, call. = FALSE)
      }
      inner <- stack[[length(stack)]]
      stack[[length(stack)]] <- NULL
      stack[[length(stack)]] <- .add_counts(stack[[length(stack)]], inner * mult)
    } else if (grepl("^[A-Z]", tk)) {
      i <- i + 1L
      k <- 1
      if (i <= n && grepl("^[0-9]", tokens[i])) {
        k <- as.numeric(tokens[i])
        i <- i + 1L
      }
      stack[[length(stack)]] <- .add_counts(stack[[length(stack)]],
                                            stats::setNames(k, tk))
    } else {
      stop("Formule non analysable : ", s, call. = FALSE)
    }
  }
  if (length(stack) != 1L) {
    stop("Parenth\u00e8ses non appari\u00e9es dans la formule : ", s, call. = FALSE)
  }
  stack[[1L]]
}

# Sorted element symbols of a parsed formula
.element_names <- function(p) {
  nm <- names(p$elements)
  if (is.null(nm)) character(0) else sort(nm)
}

# Canonical composition key, e.g. "C2H6O|0": same key = same composition and
# charge, whatever the way the formula is written.
.composition_key <- function(p) {
  nm <- names(p$elements)
  if (is.null(nm)) nm <- character(0)
  o <- order(nm)
  paste0(paste0(nm[o], p$elements[o], collapse = ""), "|", p$charge)
}

.normalize_phase <- function(phase) {
  if (is.null(phase)) return(NULL)
  map <- c(aq = "aq", aqueous = "aq", cr = "cr", solid = "cr", crystal = "cr",
           s = "cr", liq = "liq", liquid = "liq", l = "liq", gas = "gas",
           g = "gas")
  p <- map[tolower(as.character(phase))]
  if (anyNA(p)) {
    stop("Phase non reconnue. Phases accept\u00e9es : aq, cr (solide), liq, gas.",
         call. = FALSE)
  }
  unname(p)
}

#' Find the species of the database matching a query
#'
#' The query is matched, in this order, against (1) the formula written
#' exactly as in the database, (2) the name or abbreviation, ignoring case,
#' (3) the elemental composition and charge (so that `"C2H6O"` finds
#' `"C2H5OH"`). Several species may share the same formula (isomers) or the
#' same name in different phases: all of them are returned.
#'
#' @param query A single character string.
#' @param db A database table, by default [thermo_db()].
#' @return A list with `index` (row numbers of `db`) and `match` (`"formula"`,
#'   `"name"`, `"composition"` or `NA` if nothing was found).
#' @examples
#' find_species("methane")$index
#' @export
find_species <- function(query, db = thermo_db()) {
  if (!is.character(query) || length(query) != 1L || is.na(query) ||
      !nzchar(trimws(query))) {
    stop("`query` doit \u00eatre une cha\u00eene de caract\u00e8res non vide.", call. = FALSE)
  }
  q <- trimws(query)
  idx <- which(db$formula == q)
  if (length(idx)) return(list(index = idx, match = "formula"))
  ql <- tolower(q)
  idx <- which(tolower(db$name) == ql | tolower(db$abbrv) == ql)
  if (length(idx)) return(list(index = idx, match = "name"))
  key <- tryCatch(.composition_key(parse_formula(q)),
                  error = function(e) NA_character_)
  if (!is.na(key)) {
    idx <- which(db$composition == key)
    if (length(idx)) return(list(index = idx, match = "composition"))
  }
  list(index = integer(0), match = NA_character_)
}
