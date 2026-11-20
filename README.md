# standardenthalpieformationR

Package R pour obtenir l'enthalpie standard de formation ΔfH° (et les grandeurs
connexes : énergie de Gibbs de formation ΔfG°, entropie S°, capacité
calorifique Cp°, volume V°) de composés chimiques, avec leur **source**, leurs
**conditions de référence** et leur **traçabilité**.

> Version 0.1.0. Le package donne accès à la base de données et aux fonctions
> qui renvoient une information à partir d'elle. Il ne contient aucune
> méthode d'estimation : une valeur absente reste absente.

## Données

La base provient d'**OBIGT**, la base de données thermodynamiques du package R
[CHNOSZ](https://chnosz.net/) (J. Dick, licence GPL-3). 9 des 10 fichiers par
défaut de CHNOSZ sont copiés **sans modification** (3 007 espèces) dans
`inst/extdata/OBIGT/` ; `SOURCE.txt` donne la version de CHNOSZ, le commit, la
date de récupération et l'empreinte SHA-256 de chaque fichier.

| Rangement par nature (fichiers OBIGT) | Phases |
| --- | --- |
| `inorganic_aq`, `organic_aq`, `H2O_aq` | aq (espèces aqueuses), liq (eau) |
| `inorganic_cr`, `organic_cr` | cr (cristaux, polymorphes cr2 à cr9) |
| `inorganic_gas`, `organic_gas` | gas |
| `inorganic_liq`, `organic_liq` | liq |

Ce qu'il faut savoir sur les valeurs :

* La colonne `H` d'OBIGT est documentée comme « Standard molal enthalpy of
  formation from the elements » ; `G` est l'énergie de Gibbs de formation ;
  `S` et `Cp` sont des valeurs **absolues**, pas des grandeurs de formation.
* Conditions de référence : 25 °C (298,15 K) et 1 bar.
* OBIGT mélange les unités ligne par ligne (cal ou J, colonne `E_units`). Le
  package convertit tout en J/mol (1 cal = 4,184 J) et garde la valeur et
  l'unité d'origine (colonnes `orig_*`).
* OBIGT ne fournit **pas d'incertitude** : la colonne `uncertainty` est `NA`.
* Le fichier `Berman_cr` de CHNOSZ (92 minéraux sans aucune valeur G/H/S/Cp)
  n'est pas inclus : il ne contient rien à renvoyer.
* Plusieurs espèces peuvent partager une formule (isomères) : toutes sont
  renvoyées, la colonne `name` permet de choisir.

### Limites connues (constatées sur les fichiers)

* **Eau liquide** : dans OBIGT, l'entrée `water` (phase liq) n'a pas de valeur
  `H` dans le CSV (ses propriétés sont calculées par le modèle spécial de
  l'eau de CHNOSZ). `get_formation_enthalpy("H2O")` renvoie donc, pour la
  phase liquide, le statut `no_value` ; la glace (`ice`) et la vapeur
  (`steam`) ont une valeur. Cette lacune devra être comblée par une autre
  source vérifiée (CODATA, ATcT, NIST), jamais par une valeur écrite de mémoire.
* **Entrées incohérentes** : sur 2 853 entrées avec G, H et S, 66 s'écartent de
  plus de 5 kJ/mol du recalcul ΔfG° = ΔfH° − T ΔfS° (voir
  `check_gibbs_consistency()`). Elles sont signalées, pas corrigées.
* **Pas d'incertitudes** ni de SMILES/InChIKey dans OBIGT.
* **Une seule source** pour l'instant : pas de gestion de valeurs
  contradictoires entre sources tant qu'une seconde source n'est pas ajoutée
  (la détection des doublons internes, `find_duplicates()`, ne trouve rien).

## Installation

```r
# depuis le dossier du package
install.packages("devtools")          # une seule fois
devtools::document()                  # génère l'aide (man/) à partir des commentaires roxygen
devtools::install()
```

ou, depuis le zip décompressé : `devtools::install_local("standardenthalpieformationR")`.

## Utilisation

```r
library(standardenthalpieformationR)

thermo <- ThermoPredictor()

get_formation_enthalpy("methane", phase = "gas")
get_formation_enthalpy("CH4", phase = "gas", units = "kcal/mol")

# toutes les phases disponibles, sans deviner
get_formation_enthalpy("H2O")

# grandeurs connexes (ΔfH°, ΔfG°, S°, Cp°, V°)
get_thermo_properties("methane", phase = "gas")

# familles chimiques
list_families()
species_in_family("hydrocarbons", phase = "gas")

# provenance, références, contrôles
database_info()
get_references("Kel60")
check_database()
check_gibbs_consistency()
```

Chaque résultat contient : valeur, unité, phase, température et pression de
référence, méthode (toujours `database`), source, références bibliographiques,
incertitude, statut (`ok` : valeur documentée ; `no_value` : espèce présente
sans valeur ; `not_found` : espèce absente de la base) et une note.

## Familles chimiques

`list_families()` donne les 13 familles demandées. Le classement est fait par
règles à partir de la formule, de la charge, de la phase et du fichier OBIGT
(OBIGT ne contient pas les structures moléculaires). Deux familles ne sont
**pas déterminables** avec ces données et restent vides plutôt que devinées :
*organométalliques* (liaison métal-carbone) et *biomolécules* (le fichier
`biotic` a été fusionné dans `organic_aq`). Les radicaux ne sont pas
identifiables non plus : la famille « Ions et radicaux » ne contient que des
ions.

## Contrôle de cohérence des données

`check_gibbs_consistency()` recalcule ΔfG° = ΔfH° − T (S − Σ S°éléments) pour
chaque entrée qui possède G, H et S, avec les entropies des éléments fournies
par CHNOSZ. Les écarts supérieurs à la tolérance (5 kJ/mol par défaut) sont
signalés : ils désignent des entrées à examiner, pas forcément des erreurs.

## Organisation du code

```
R/
  chemical_identification.R   analyse des formules, recherche d'espèces
  import_data.R               lecture et normalisation d'OBIGT
  families.R                  familles chimiques
  database.R                  accès à la base, références, doublons
  prediction.R                résultat, propriétés connexes
  method_selection.R          ThermoPredictor, get_formation_enthalpy()
  validation.R                contrôles d'intégrité et de cohérence
  units.R                     conversions d'unités
```

Une valeur renvoyée avec `method_class = "database"` est une valeur
documentée de la base. Le package ne produit aucune estimation : une valeur
absente de la base reste absente (statuts `no_value` et `not_found`).

## Vérifier le package

R n'était pas disponible dans l'environnement où ce package a été écrit : le
code n'a **pas encore été exécuté**. Sa logique d'analyse des formules et les
comptages des tests ont été vérifiés avec un programme équivalent sur les
3 007 entrées. Pour valider chez vous :

```r
devtools::document()
devtools::test()
devtools::check()
```

Si une erreur apparaît, copiez le message : il suffit pour corriger.

## Mettre à jour les données

Voir `data-raw/update_obigt.R` (copie des fichiers depuis une copie locale de
CHNOSZ), puis mettre à jour `inst/extdata/OBIGT/SOURCE.txt`.

## Licence et crédits

GPL-3 (comme CHNOSZ, dont les données sont reprises). Si vous utilisez ces
données, citez CHNOSZ : Dick, J.M. (2019) *CHNOSZ: Thermodynamic calculations
and diagrams for geochemistry*, Frontiers in Earth Science 7:180,
doi:10.3389/feart.2019.00180. Les références de chaque valeur sont dans les
colonnes `ref1` et `ref2` (voir `get_references()`).

Dans `DESCRIPTION`, remplacez l'adresse `REPLACE_ME@example.com` par la vôtre.
