#' @title Resumen de Modelo RegenerAR
#' @description Genera un resumen estadístico para modelos de RegenerAR.
#' @param object Objeto de clase modelo_germinacion_unificado.
#' @param ... Argumentos adicionales.
#' @return Objeto de clase summary.modelo_germinacion_unificado.
#' @export
summary.modelo_germinacion_unificado <- function(object, ...) {
  modelo_base <- object$modelo
  tipo <- object$tipo_modelo
  familia <- object$familia

  if (inherits(modelo_base, "glm") || inherits(modelo_base, "lm")) {
    res_base <- summary(modelo_base)
  } else {
    res_base <- list(coefficients = coef(modelo_base), call = object$call)
  }

  out <- list(
    call = res_base$call,
    familia = familia,
    tipo = tipo,
    formula = object$formula,
    coefficients = res_base$coefficients,
    aic = if (!is.null(modelo_base$aic)) modelo_base$aic else NA,
    deviance = if (!is.null(modelo_base$deviance)) modelo_base$deviance else NA,
    null.deviance = if (!is.null(modelo_base$null.deviance)) modelo_base$null.deviance else NA,
    df.residual = modelo_base$df.residual,
    diagnostico = object$diagnostico,
    factores_activos = object$factores_activos
  )

  class(out) <- "summary.modelo_germinacion_unificado"
  return(out)
}

#' @title Imprimir Resumen de Modelo
#' @description Imprime el resumen de un modelo RegenerAR.
#' @param x Objeto de clase summary.modelo_germinacion_unificado.
#' @param digits Número de decimales.
#' @param signif.stars Lógico. Mostrar estrellas de significancia.
#' @param ... Argumentos adicionales.
#' @export
print.summary.modelo_germinacion_unificado <- function(x, digits = 4, signif.stars = TRUE, ...) {

  cat("\n=== Resumen del Modelo RegenerAR ===\n\n")
  cat("Fórmula:", deparse(x$formula), "\n")
  cat("Familia:", x$familia, "| Estructura:", x$tipo, "\n")
  cat("Factores:", paste(x$factores_activos, collapse = ", "), "\n\n")

  coefs <- x$coefficients

  if (ncol(coefs) >= 4) {
    p_vals <- coefs[, 4]

    # Calcular estrellas de significancia
    stars <- if (signif.stars) {
      symnum(p_vals, corr = FALSE, na = FALSE,
             cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
             symbols = c("***", "**", "*", ".", " "))
    } else {
      rep("", length(p_vals))
    }

    # Crear la tabla de coeficientes
    tabla_print <- data.frame(
      Estimate = formatC(coefs[, 1], format = "f", digits = digits),
      `Std. Error` = formatC(coefs[, 2], format = "f", digits = digits),
      `z/t value` = formatC(coefs[, 3], format = "f", digits = digits),
      `Pr(>|z|)` = format.pval(p_vals, digits = digits),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    # 🔧 FIX: Asignar estrellas de forma segura (sin usar $``)
    if (signif.stars) {
      tabla_print$stars <- stars
      names(tabla_print)[ncol(tabla_print)] <- "" # Renombrar la última columna a vacío
    }

    print(tabla_print, right = TRUE)

    if (signif.stars) {
      cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
    }
  } else {
    print(coefs)
  }

  # Diagnósticos y Ajuste
  if (!is.na(x$aic)) cat("AIC:", round(x$aic, 2), "\n")
  if (!is.na(x$deviance)) {
    cat("Deviance Residuals (Residual):", round(x$deviance, 2),
        "on", x$df.residual, "degrees of freedom\n")
  }
  if (!is.na(x$null.deviance)) {
    cat("Null Deviance:", round(x$null.deviance, 2), "\n")
  }

  # Mostrar diagnósticos de RegenerAR si existen
  if (length(x$diagnostico) > 0) {
    cat("\n[Diagnósticos RegenerAR]\n")
    if (!is.null(x$diagnostico$dispersion)) {
      cat("- Dispersión: Ratio =", round(x$diagnostico$dispersion$ratio, 2),
          "(p =", format.pval(x$diagnostico$dispersion$p_value), ")\n")
    }
    if (!is.null(x$diagnostico$separacion) && isTRUE(x$diagnostico$separacion$detectada)) {
      cat("- Separación perfecta detectada en:\n")
      for (comb in x$diagnostico$separacion$combinaciones) {
        cat("   *", comb, "\n")
      }
    }
  }

  invisible(x)
}

