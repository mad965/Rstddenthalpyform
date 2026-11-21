# R/prediction_model.R
# Modèles de régression multiple pour ΔfH° (un par source, sans filtre de phase).
# ΔfH° = Σ β_X · n_X

# ---------------------------------------------------------------------------
# 1. Utilitaires
# ---------------------------------------------------------------------------

#' Compter les occurrences de chaque élément dans les formules
#'
#' @param df data.frame avec une colonne `formula`.
#' @return data.frame `element`, `n_composes`, trié par fréquence.
#' @keywords internal
count_elements <- function(df) {
  elems <- lapply(df$formula, function(f) {
    tryCatch(names(parse_formula(f)$elements), error = function(e) character(0))
  })
  all_e <- unlist(elems)
  res <- as.data.frame(table(all_e), stringsAsFactors = FALSE)
  names(res) <- c("element", "n_composes")
  res[order(-res$n_composes), , drop = FALSE]
}

#' Construire la matrice X des counts d'éléments
#'
#' @param df data.frame avec `formula`, `delta_fH`.
#' @param elements Vecteur de noms d'éléments (colonnes de X).
#' @return Liste avec X, y, df (lignes conservées).
#' @keywords internal
build_design_matrix <- function(df, elements) {
  ok <- !is.na(df$delta_fH) & !is.na(df$formula) & nzchar(df$formula)
  df <- df[ok, , drop = FALSE]
  
  parsed <- lapply(df$formula, function(f) {
    tryCatch(parse_formula(f)$elements, error = function(e) NULL)
  })
  keep <- !vapply(parsed, is.null, logical(1))
  df <- df[keep, , drop = FALSE]
  parsed <- parsed[keep]
  
  X <- matrix(0, nrow = length(parsed), ncol = length(elements),
              dimnames = list(NULL, elements))
  for (i in seq_along(parsed)) {
    cnt <- parsed[[i]]
    common <- intersect(names(cnt), elements)
    if (length(common) > 0) X[i, common] <- cnt[common]
  }
  
  list(X = X, y = df$delta_fH, df = df)
}

# ---------------------------------------------------------------------------
# 2. Ajustement
# ---------------------------------------------------------------------------

#' Ajuster un modèle linéaire de ΔfH° sur les counts d'éléments
#'
#' @param df data.frame avec `formula`, `delta_fH`.
#' @param elements Vecteur d'éléments (prédicteurs).
#' @param intercept Inclure un intercept ? Défaut FALSE.
#' @return Un objet lm avec attributs `elements` et `n_train`.
#' @export
fit_enthalpy_model <- function(df, elements, intercept = FALSE) {
  dm <- build_design_matrix(df, elements)
  
  if (nrow(dm$X) < length(elements) + 1) {
    stop("Pas assez de composés (", nrow(dm$X), ") pour ",
         length(elements), " prédicteurs.", call. = FALSE)
  }
  
  data_fit <- as.data.frame(dm$X)
  data_fit$y <- dm$y
  
  form <- if (intercept) {
    as.formula(paste("y ~", paste(elements, collapse = " + ")))
  } else {
    as.formula(paste("y ~ 0 +", paste(elements, collapse = " + ")))
  }
  
  model <- lm(form, data = data_fit)
  attr(model, "elements") <- elements
  attr(model, "n_train")  <- nrow(dm$X)
  model
}

# ---------------------------------------------------------------------------
# 3. Prédiction
# ---------------------------------------------------------------------------

#' Prédire ΔfH° pour un composé à partir de sa formule
#'
#' @param model Objet lm retourné par fit_enthalpy_model().
#' @param formula Formule brute (ex. "CH4").
#' @return Valeur prédite en kJ/mol (numérique).
#' @export
predict_enthalpy <- function(model, formula) {
  elements <- attr(model, "elements")
  cnt <- tryCatch(parse_formula(formula)$elements, error = function(e) NULL)
  if (is.null(cnt)) return(NA_real_)
  newdata <- as.data.frame(as.list(setNames(rep(0, length(elements)), elements)))
  common <- intersect(names(cnt), elements)
  if (length(common) > 0) newdata[1, common] <- cnt[common]
  as.numeric(predict(model, newdata = newdata))
}

# ---------------------------------------------------------------------------
# 4. Ajustement des 3 modèles
# ---------------------------------------------------------------------------

#' Ajuster les modèles ATcT, RMG et CHNOSZ avec les listes d'éléments retenues
#'
#' @param atct,rmg,chnosz data.frames harmonisés.
#' @return Liste nommée de 3 objets lm.
#' @export
fit_all_models <- function(atct, rmg, chnosz) {
  list(
    ATcT = fit_enthalpy_model(atct,
                              elements = c("H","C","O","N","Cl","F","Br","I","Si","S")),
    RMG = fit_enthalpy_model(rmg,
                             elements = c("C","H","O","F","Br","Cl","N","S","Li")),
    CHNOSZ = fit_enthalpy_model(chnosz,
                                elements = c("H","C","O","N","S","P"))
  )
}

# ---------------------------------------------------------------------------
# 5. Comparaison des coefficients
# ---------------------------------------------------------------------------

#' Tableau comparatif des coefficients des 3 modèles
#'
#' @param models Liste retournée par fit_all_models().
#' @return data.frame avec une ligne par élément, une colonne par source.
#' @export
compare_model_coefficients <- function(models) {
  all_e <- sort(unique(unlist(lapply(models, function(m) attr(m, "elements")))))
  out <- data.frame(element = all_e, stringsAsFactors = FALSE)
  for (src in names(models)) {
    out[[src]] <- coef(models[[src]])[out$element]
  }
  out
}

# ---------------------------------------------------------------------------
# 6. Graphique prédit vs réel (avec couleur par ligne du tableau périodique
#    et dégradé par état physique)
# ---------------------------------------------------------------------------

#' Numéro de ligne (période) d'un élément dans le tableau périodique
#'
#' @param element Symbole chimique (ex. "C").
#' @return Entier 1–7, ou NA si inconnu.
#' @keywords internal
element_period <- function(element) {
  periods <- list(
    `1` = c("H","He"),
    `2` = c("Li","Be","B","C","N","O","F","Ne"),
    `3` = c("Na","Mg","Al","Si","P","S","Cl","Ar"),
    `4` = c("K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
            "Ga","Ge","As","Se","Br","Kr"),
    `5` = c("Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
            "In","Sn","Sb","Te","I","Xe"),
    `6` = c("Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy",
            "Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt",
            "Au","Hg","Tl","Pb","Bi","Po","At","Rn"),
    `7` = c("Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf",
            "Es","Fm","Md","No","Lr")
  )
  for (p in names(periods)) {
    if (element %in% periods[[p]]) return(as.integer(p))
  }
  NA_integer_
}

#' Graphique prédit vs réel pour un modèle
#'
#' Trace les points (ΔfH° réel, ΔfH° prédit), colorés par ligne du tableau
#' périodique (teinte) et par état physique (dégradé de la même teinte).
#' Ajoute la droite y = x.
#'
#' @param model Objet lm retourné par fit_enthalpy_model().
#' @param df data.frame utilisé pour l'ajustement (ou un jeu de test).
#' @param titre Titre du graphique.
#' @return Un objet ggplot (invisiblement).
#' @export
plot_predicted_vs_true <- function(model, df, titre = "Prédit vs réel") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Le package ggplot2 est requis.", call. = FALSE)
  }
  elements <- attr(model, "elements")
  
  # Prédictions
  dm <- build_design_matrix(df, elements)
  pred <- as.numeric(predict(model, newdata = as.data.frame(dm$X)))
  df_ok <- dm$df
  df_ok$pred <- pred
  
  # Ligne du tableau périodique : on prend le premier élément de la formule
  # qui appartient à `elements` (le plus "lourd" serait plus juste, mais
  # pour un premier jet, on prend le premier trouvé).
  df_ok$period <- vapply(df_ok$formula, function(f) {
    cnt <- tryCatch(parse_formula(f)$elements, error = function(e) NULL)
    if (is.null(cnt)) return(NA_integer_)
    common <- intersect(names(cnt), elements)
    if (length(common) == 0) return(NA_integer_)
    element_period(common[1])
  }, integer(1))
  
  df_ok$period_f <- factor(df_ok$period, levels = 1:7,
                           labels = paste0("Période ", 1:7))
  
  ggplot2::ggplot(df_ok, ggplot2::aes(x = delta_fH, y = pred,
                                      colour = period_f)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey40") +
    ggplot2::geom_point(ggplot2::aes(alpha = phase), size = 2) +
    ggplot2::scale_colour_brewer(palette = "Set1", na.value = "grey70") +
    ggplot2::scale_alpha_manual(values = c(g = 0.9, l = 0.6, cr = 0.4,
                                           aq = 0.3), na.value = 0.5) +
    ggplot2::labs(
      title = titre,
      x = expression(Delta[f]*H*"° réel (kJ/mol)"),
      y = expression(Delta[f]*H*"° prédit (kJ/mol)"),
      colour = "Ligne (TCPE)",
      alpha  = "Phase"
    ) +
    ggplot2::theme_minimal()
}