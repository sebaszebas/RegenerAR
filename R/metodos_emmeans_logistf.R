# ==============================================================================
# MÉTODOS S3 PARA COMPATIBILIDAD logistf <-> emmeans
# Permiten que emmeans funcione con modelos ajustados por logistf (Firth)
# ==============================================================================

#' Recuperar datos de un modelo logistf para emmeans
#'
#' @param object Objeto de clase \code{logistf}
#' @param trms Términos del modelo
#' @param xlev Niveles de los factores
#' @param grid Grid de referencia (no usado aquí)
#' @param ... Argumentos adicionales
#'
#' @return Objeto de datos recuperados compatible con emmeans
#' @keywords internal
#' @exportS3Method emmeans::recover_data
recover_data.logistf <- function(object, trms, xlev, grid, ...) {
  fcall <- object$call

  # Recuperar el entorno real de la fórmula
  env <- environment(object$formula)
  if (is.null(env)) env <- parent.frame()

  # Evaluar el data frame del modelo
  data <- tryCatch(
    eval(fcall$data, envir = env),
    error = function(e) NULL
  )

  # Fallback: si es símbolo o carácter, buscarlo con get()
  if (is.null(data) || is.symbol(fcall$data) || is.character(fcall$data)) {
    data_name <- as.character(fcall$data)
    data <- tryCatch(
      get(data_name, envir = env),
      error = function(e) {
        cli::cli_warn(c(
          "!" = "No se pudo recuperar los datos del modelo logistf.",
          "i" = "emmeans podría no funcionar correctamente."
        ))
        return(NULL)
      }
    )
  }

  if (is.null(data)) {
    cli::cli_abort("No se pudieron recuperar los datos para el modelo logistf.")
  }

  # Llamar al método genérico de emmeans
  emmeans::recover_data(fcall, trms, xlev, data = data, ...)
}

#' Base para medias marginales estimadas de un modelo logistf
#'
#' @param object Objeto de clase \code{logistf}
#' @param trms Términos del modelo
#' @param xlev Niveles de los factores
#' @param grid Grid de referencia donde evaluar las medias
#' @param ... Argumentos adicionales
#'
#' @return Lista con la base para calcular emmeans
#' @keywords internal
#' @exportS3Method emmeans::emm_basis
emm_basis.logistf <- function(object, trms, xlev, grid, ...) {
  # Extraer coeficientes y matriz de varianza-covarianza
  b <- stats::coef(object)
  V <- stats::vcov(object)

  # Recuperar datos usando nuestro método blindado
  data <- recover_data.logistf(object, trms, xlev, grid)

  # Construir matriz de diseño manualmente (logistf no tiene model.matrix nativo)
  X <- stats::model.matrix(object$formula, data = data)

  #  CRÍTICO: Alinear nombres para evitar desajustes numéricos
  common_names <- intersect(names(b), colnames(X))
  if (length(common_names) == 0) {
    cli::cli_abort(c(
      "x" = "No hay coincidencia entre coeficientes y matriz de diseño.",
      "i" = "Coeficientes: {.val {names(b)}}",
      "i" = "Columnas X: {.val {colnames(X)}}"
    ))
  }

  # Filtrar a nombres comunes
  b <- b[common_names]
  X <- X[, common_names, drop = FALSE]
  V <- V[common_names, common_names, drop = FALSE]

  # Construir la matriz de diseño para el grid de emmeans
  X_grid <- stats::model.matrix(trms, data = grid, xlev = xlev)

  # Alinear con los coeficientes
  X_grid <- X_grid[, common_names, drop = FALSE]

  # Para Firth, usamos grados de libertad infinitos (aproximación Wald)
  list(
    X = X_grid,
    b = b,
    V = V,
    dffun = function(k, dfargs) Inf,
    dfargs = list(),
    misc = list(
      link = "logit",
      fam = "binomial"
    )
  )
}
