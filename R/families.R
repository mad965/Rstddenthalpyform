# Chemical families. The classification is rule-based and derived only from
# the formula, the charge, the phase and the OBIGT file of each entry:
# OBIGT does not store molecular structures. Families that cannot be
# determined from these data are declared with `determinable = FALSE` and
# left empty rather than guessed.
.family_defs <- data.frame(
  family = c("elements", "inorganic", "organic", "hydrocarbons", "oxygenated",
             "nitrogenous", "sulfur", "halogenated", "organometallic",
             "minerals", "aqueous", "ions", "biomolecules"),
  label = c("\u00c9l\u00e9ments", "Compos\u00e9s inorganiques", "Compos\u00e9s organiques",
            "Hydrocarbures", "Compos\u00e9s oxyg\u00e9n\u00e9s", "Compos\u00e9s azot\u00e9s",
            "Compos\u00e9s soufr\u00e9s", "Compos\u00e9s halog\u00e9n\u00e9s",
            "Compos\u00e9s organom\u00e9talliques", "Min\u00e9raux", "Esp\u00e8ces aqueuses",
            "Ions et radicaux", "Biomol\u00e9cules"),
  rule = c(
    "Formule d'un seul \u00e9l\u00e9ment, neutre, phase non aqueuse (inclut les allotropes hors \u00e9tat de r\u00e9f\u00e9rence)",
    "Classement d'OBIGT : fichiers autres que organic_*",
    "Classement d'OBIGT : fichiers organic_*",
    "Formule limit\u00e9e \u00e0 C et H, neutre",
    "Organique contenant C et O",
    "Organique contenant C et N",
    "Organique contenant C et S",
    "Organique contenant C et F, Cl, Br ou I",
    "Non d\u00e9terminable : OBIGT ne contient pas la structure (liaison m\u00e9tal-carbone)",
    "Phase cristalline inorganique de plus d'un \u00e9l\u00e9ment",
    "Phase aqueuse (aq)",
    "Charge non nulle ; les radicaux ne sont pas identifiables dans OBIGT",
    "Non d\u00e9terminable : le fichier biotic a \u00e9t\u00e9 fusionn\u00e9 dans organic_aq (CHNOSZ 2.0.0)"
  ),
  determinable = c(rep(TRUE, 8), FALSE, rep(TRUE, 3), FALSE),
  stringsAsFactors = FALSE
)

.build_family_table <- function(db) {
  els <- strsplit(db$elements, " ", fixed = TRUE)
  n_el <- lengths(els)
  has <- function(sym) vapply(els, function(e) sym %in% e, logical(1))
  only <- function(set) {
    vapply(els, function(e) length(e) > 0L && all(e %in% set), logical(1))
  }
  halo <- vapply(els, function(e) any(e %in% c("F", "Cl", "Br", "I")),
                 logical(1))
  nat <- db$nature
  ph <- db$phase
  chg <- db$charge
  none <- rep(FALSE, nrow(db))
  fam <- list(
    elements = n_el == 1L & chg == 0 & ph != "aq",
    inorganic = nat == "inorganic",
    organic = nat == "organic",
    hydrocarbons = has("C") & has("H") & only(c("C", "H")) & chg == 0,
    oxygenated = nat == "organic" & has("C") & has("O"),
    nitrogenous = nat == "organic" & has("C") & has("N"),
    sulfur = nat == "organic" & has("C") & has("S"),
    halogenated = nat == "organic" & has("C") & halo,
    organometallic = none,
    minerals = ph == "cr" & nat == "inorganic" & n_el > 1L,
    aqueous = ph == "aq",
    ions = chg != 0,
    biomolecules = none
  )
  rows <- lapply(names(fam), function(f) {
    i <- which(fam[[f]])
    if (length(i)) {
      data.frame(id = db$id[i], family = rep(f, length(i)),
                 stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}
