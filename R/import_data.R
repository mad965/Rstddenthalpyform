# OBIGT files loaded by default by CHNOSZ (same selection as the "default"
# files of the CHNOSZ vignette on OBIGT). Optional files (SUPCRT92, SLOP98,
# DEW, AD) are not included: they partly duplicate the default species with
# other values. Berman_cr is not included either: its 92 minerals have no
# G, H, S or Cp value in OBIGT.
.obigt_files <- c(
  "H2O_aq", "inorganic_aq", "organic_aq", "inorganic_cr", "organic_cr",
  "inorganic_gas", "organic_gas", "inorganic_liq", "organic_liq"
)

.required_cols <- c("name", "abbrv", "formula", "state", "ref1", "ref2",
                    "date", "model", "E_units", "G", "H", "S", "Cp", "V")

.obigt_dir <- function() {
  p <- system.file("extdata", "OBIGT", package = "Rstdenthalpyform")
  if (nzchar(p) && dir.exists(p)) return(p)
  if (dir.exists("inst/extdata/OBIGT")) return("inst/extdata/OBIGT")
  ""
}

.read_source_info <- function(dir = .obigt_dir()) {
  path <- file.path(dir, "SOURCE.txt")
  if (!nzchar(dir) || !file.exists(path)) return(list())
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#")]
  key <- sub(":.*$", "", lines)
  val <- trimws(sub("^[^:]*:", "", lines))
  stats::setNames(as.list(val), key)
}

.source_label <- function(dir = .obigt_dir()) {
  info <- .read_source_info(dir)
  ver <- info[["chnosz_version"]]
  if (is.null(ver)) return("OBIGT (CHNOSZ)")
  paste0("OBIGT (CHNOSZ ", sub(" .*$", "", ver), ")")
}

#' Import the OBIGT files into a normalised table
#'
#' Reads the OBIGT CSV files shipped in `inst/extdata/OBIGT` and builds one
#' normalised table. Nothing is modified in the raw files: every value is
#' kept in its original unit (`orig_*` columns) next to its conversion to
#' joules (thermochemical calorie, 1 cal = 4.184 J). A missing value stays
#' `NA` and is never replaced by zero.
#'
#' Columns of the result: `id`, `name`, `abbrv`, `formula`, `nature`
#' (`"inorganic"` or `"organic"`, from the name of the OBIGT file), `phase`
#' (`"aq"`, `"cr"`, `"liq"` or `"gas"`), `state` (original OBIGT state, e.g.
#' `"cr2"` for a second polymorph), `dfG_J_mol`, `dfH_J_mol`, `S_J_mol_K`,
#' `Cp_J_mol_K`, `V_cm3_mol`, the `orig_*` columns, `T_ref_K`, `p_ref_bar`,
#' `model`, `ref1`, `ref2`, `entry_date`, `source`, `source_file`,
#' `source_row`, `charge`, `elements` and `composition`.
#'
#' @param dir Directory containing the OBIGT CSV files.
#' @return A data frame with one row per OBIGT entry.
#' @export
import_obigt <- function(dir = .obigt_dir()) {
  if (!nzchar(dir) || !dir.exists(dir)) {
    stop("Dossier de donn\u00e9es OBIGT introuvable : '", dir, "'.", call. = FALSE)
  }
  parts <- lapply(.obigt_files, function(stem) {
    path <- file.path(dir, paste0(stem, ".csv"))
    if (!file.exists(path)) {
      stop("Fichier OBIGT manquant : ", path, call. = FALSE)
    }
    raw <- utils::read.csv(path, colClasses = "character",
                           na.strings = c("NA", ""), check.names = FALSE,
                           encoding = "UTF-8", stringsAsFactors = FALSE)
    absent <- setdiff(.required_cols, names(raw))
    if (length(absent)) {
      stop("Colonnes absentes de ", basename(path), " : ",
           paste(absent, collapse = ", "), call. = FALSE)
    }
    raw <- raw[, .required_cols, drop = FALSE]
    raw$source_file <- rep(stem, nrow(raw))
    raw$source_row <- seq_len(nrow(raw))
    raw
  })
  raw <- do.call(rbind, parts)
  rownames(raw) <- NULL

  eu <- raw$E_units
  if (anyNA(eu) || !all(eu %in% c("J", "cal"))) {
    stop("Unit\u00e9s d'\u00e9nergie inattendues dans OBIGT (attendu : J ou cal).",
         call. = FALSE)
  }
  if (anyNA(raw$formula) || anyNA(raw$name) || anyNA(raw$state)) {
    stop("Entr\u00e9es OBIGT sans nom, formule ou \u00e9tat.", call. = FALSE)
  }
  fac <- ifelse(eu == "cal", .CAL_J, 1)
  num <- function(x) suppressWarnings(as.numeric(x))
  g <- num(raw$G)
  h <- num(raw$H)
  s <- num(raw$S)
  cp <- num(raw$Cp)

  phase <- sub("[0-9]+$", "", raw$state)
  if (!all(phase %in% c("aq", "cr", "liq", "gas"))) {
    stop("\u00c9tats OBIGT inattendus.", call. = FALSE)
  }

  parsed <- lapply(raw$formula, parse_formula)

  out <- data.frame(
    id = sprintf("OBIGT:%s:%d", raw$source_file, raw$source_row),
    name = raw$name,
    abbrv = raw$abbrv,
    formula = raw$formula,
    nature = ifelse(startsWith(raw$source_file, "organic"),
                    "organic", "inorganic"),
    phase = phase,
    state = raw$state,
    dfG_J_mol = g * fac,
    dfH_J_mol = h * fac,
    S_J_mol_K = s * fac,
    Cp_J_mol_K = cp * fac,
    V_cm3_mol = num(raw$V),
    orig_unit = eu,
    orig_dfG = g,
    orig_dfH = h,
    orig_S = s,
    orig_Cp = cp,
    T_ref_K = 298.15,
    p_ref_bar = 1,
    model = raw$model,
    ref1 = raw$ref1,
    ref2 = raw$ref2,
    entry_date = raw$date,
    source = .source_label(dir),
    source_file = raw$source_file,
    source_row = raw$source_row,
    charge = vapply(parsed, function(p) p$charge, numeric(1)),
    elements = vapply(parsed, function(p) {
      paste(.element_names(p), collapse = " ")
    }, character(1)),
    composition = vapply(parsed, .composition_key, character(1)),
    stringsAsFactors = FALSE
  )
  out
}
