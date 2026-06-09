##############################################################################
# splitBy_local.R — Reemplazo local de doBy::splitBy()
#
# Elimina la dependencia de doBy (que arrastra Deriv, broom, modelr,
# microbenchmark, etc. y genera conflictos de versiones con r-base/ggplot2).
#
# Solo implementa el caso que usan los scripts de sppIDer:
#   splitBy("column_name", dataframe)
#
# Devuelve una lista nombrada de dataframes, igual que doBy::splitBy.
# Licencia: GPL-2 (compatible con doBy original).
##############################################################################

splitBy <- function(formula, data = parent.frame(), drop = TRUE) {
  # Caso 1: formula es un string -> nombre de columna directamente
  if (is.character(formula)) {
    vars <- formula
  }
  # Caso 2: formula es una formula como ~species
  else if (inherits(formula, "formula")) {
    vars <- all.vars(formula)
  }
  else {
    stop("splitBy: 'formula' debe ser un string o una formula (~var)")
  }

  if (length(vars) != 1) {
    stop("splitBy local: solo soporta split por UNA variable, recibidas: ",
         paste(vars, collapse = ", "))
  }

  col <- vars[1]
  if (!col %in% names(data)) {
    stop("splitBy: columna '", col, "' no encontrada en el dataframe")
  }

  # Convertir a factor para mantener el orden de aparicion
  f <- factor(data[[col]], levels = unique(data[[col]]))

  # split() de base R
  result <- split(data, f, drop = drop)

  result
}
