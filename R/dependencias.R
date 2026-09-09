#' Instalar todas las dependencias opcionales de RegenerAR
#'
#' @description
#' Instala los paquetes sugeridos que habilitan funcionalidades avanzadas:
#' \itemize{
#'   \item \code{logistf}: Regresión logística de Firth (para separación perfecta)
#'   \item \code{performance}: Diagnóstico de sobredispersión en GLM binomial
#'   \item \code{broom}: Limpieza de outputs de modelos para tablas
#' }
#'
#' @details
#' Esta función es opcional. Si no la ejecutas, el paquete funcionará igual,
#' pero algunas funcionalidades avanzadas estarán deshabilitadas.
#'
#' @param repos Vector de caracteres con los repositorios CRAN. Por defecto
#'   usa el repositorio oficial de R.
#'
#' @return Invisiblemente, un vector lógico indicando qué paquetes fueron
#'   instalados exitosamente.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' instalar_dependencias_opcionales()
#' }
instalar_dependencias_opcionales <- function(
    repos = "https://cran.r-project.org"
) {
  paquetes <- c("logistf", "performance", "broom")

  cli::cli_inform(c(
    "i" = "Instalando dependencias opcionales de RegenerAR...",
    "i" = "Paquetes: {.val {paquetes}}"
  ))

  resultados <- stats::setNames(logical(length(paquetes)), paquetes)

  for (pkg in paquetes) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_inform("i" = "{.pkg {pkg}} ya está instalado. Saltando.")
      resultados[pkg] <- TRUE
      next
    }

    tryCatch({
      utils::install.packages(pkg, repos = repos)

      if (requireNamespace(pkg, quietly = TRUE)) {
        cli::cli_inform("v" = "{.pkg {pkg}} instalado exitosamente.")
        resultados[pkg] <- TRUE
      } else {
        cli::cli_warn("x" = "{.pkg {pkg}} no se pudo cargar tras la instalación.")
        resultados[pkg] <- FALSE
      }
    }, error = function(e) {
      cli::cli_warn(c(
        "x" = "No se pudo instalar {.pkg {pkg}}: {e$message}",
        "i" = "Intenta instalarlo manualmente con: {.code install.packages('{pkg}')}"
      ))
      resultados[pkg] <- FALSE
    })
  }

  cli::cli_inform(c(
    "v" = "Proceso completado.",
    "i" = "Reinicia R si es necesario para que los cambios surtan efecto."
  ))

  invisible(resultados)
}
