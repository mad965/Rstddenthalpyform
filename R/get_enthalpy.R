# R/get_enthalpy.R

#' Récupérer les propriétés thermodynamiques d'un composé dans les 3 sources
#'
#' Interroge ATcT, RMG et CHNOSZ (séparément, sans fusion) et retourne
#' une ligne par observation trouvée, avec delta_fH, delta_fG, S, Cp et
#' les métadonnées (phase, température, incertitude, référence).
#'
#' @param query Nom, formule ou CAS du composé.
#' @param source Vecteur de sources à interroger : "ATcT", "RMG", "CHNOSZ"
#'   ou "all".
#' @param phase Filtrer sur une phase ("g", "l", "cr", "aq"). NULL = toutes.
#' @param exact Si TRUE (défaut), correspondance exacte sur le nom ou la
#'   formule, insensible à la casse. Si FALSE, correspondance partielle
#'   (regex insensible à la casse).
#' @return Un data.frame avec une ligne par (source, phase, température)
#'   trouvée. Colonnes : source, name, formula, CAS, phase,
#'   compound_nature, delta_fH, delta_fG, S, Cp, temperature,
#'   uncertainty, unit, source_id, reference.
#' @export
get_enthalpy <- function(query,
                         source = c("ATcT", "RMG", "CHNOSZ"),
                         phase  = NULL,
                         exact  = TRUE) {
  
  source <- match.arg(source, several.ok = TRUE)
  if ("all" %in% source) source <- c("ATcT", "RMG", "CHNOSZ")
  
  # Charge un .rda à la demande (package installé OU data/ local)
  .load_one <- function(s) {
    obj_name <- paste0(tolower(s), "_harmonized")
    e <- new.env(parent = emptyenv())
    
    # 1. Chercher dans le package installé
    pkg_rda <- system.file("data", paste0(s, ".rda"),
                           package = "Rstdenthalpyform")
    if (nzchar(pkg_rda) && file.exists(pkg_rda)) {
      load(pkg_rda, envir = e)
      if (exists(obj_name, envir = e)) return(get(obj_name, envir = e))
    }
    
    # 2. Fallback : fichier local dans data/ (mode dev)
    local_rda <- file.path("data", paste0(s, ".rda"))
    if (file.exists(local_rda)) {
      load(local_rda, envir = e)
      if (exists(obj_name, envir = e)) return(get(obj_name, envir = e))
    }
    
    stop("Impossible de charger ", s, ".rda ",
         "(ni dans le package installé, ni dans data/).", call. = FALSE)
  }
  
  # Prédicat de correspondance
  .match <- function(df, q) {
    if (exact) {
      q_low <- tolower(q)
      hit <- (!is.na(df$name)    & tolower(df$name)    == q_low) |
        (!is.na(df$formula) & tolower(df$formula) == q_low)
      if (!is.null(df$CAS)) {
        hit <- hit | (!is.na(df$CAS) & df$CAS == q)
      }
      hit
    } else {
      pat <- paste0("(?i)", q)
      hit <- (!is.na(df$name)    & grepl(pat, df$name,    perl = TRUE)) |
        (!is.na(df$formula) & grepl(pat, df$formula, perl = TRUE))
      if (!is.null(df$CAS)) {
        hit <- hit | (!is.na(df$CAS) & grepl(pat, df$CAS, perl = TRUE))
      }
      hit
    }
  }
  
  # Colonnes de sortie (ordre fixe)
  cols_out <- c(
    "source", "name", "formula", "CAS", "phase", "compound_nature",
    "delta_fH", "delta_fG", "S", "Cp",
    "temperature", "uncertainty", "unit",
    "source_id", "reference"
  )
  
  # Interroge chaque source séparément
  res_list <- lapply(source, function(s) {
    df <- .load_one(s)
    hits <- df[.match(df, query), , drop = FALSE]
    if (nrow(hits) == 0) return(NULL)
    hits$source <- s
    hits[, cols_out, drop = FALSE]
  })
  
  res <- do.call(rbind, res_list)
  if (is.null(res) || nrow(res) == 0) {
    message("Aucun composé trouvé pour : ", query)
    return(res)
  }
  
  # Filtre phase
  if (!is.null(phase)) {
    res <- res[!is.na(res$phase) & res$phase %in% phase, , drop = FALSE]
    if (nrow(res) == 0) {
      message("Aucune correspondance pour la phase : ",
              paste(phase, collapse = ", "))
      return(res)
    }
  }
  
  rownames(res) <- NULL
  res
}