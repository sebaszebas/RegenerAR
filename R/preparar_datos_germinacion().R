#' Preparación Integrada de Métricas de Germinación (Matriz Maestra)
#'
#' Recibe la matriz base en formato ancho (con columnas tipo D1, D2, D3...
#' representando conteos diarios de semillas germinadas) y calcula un conjunto
#' completo de métricas de germinación: porcentaje de germinación (PG), tiempos
#' percentiles (t0, t25, t50, t75), e índices cinéticos clásicos (TMG, IVG, CVG,
#' UNC, SYN, VAR_TG, VGE, Z_INDEX, GRI), además del índice poblacional RLG.
#'
#' @details
#' ## Blindaje de datos
#' La función incluye validaciones estrictas contra valores indeterminados
#' (NaN, Inf), conversión silenciosa de tipos, y detección automatizada de
#' patrones sospechosos (datos acumulados en lugar de conteos diarios).
#'
#' ## Freno biológico
#' Si \code{solo_luz_cinetica = TRUE}, los índices cinéticos y percentiles
#' se calculan únicamente para la condición de Luz, dejando \code{NA} en
#' las filas de Oscuridad. Esto es útil cuando la oscuridad se evaluó con
#' baja frecuencia de muestreo y los índices cinéticos no son confiables.
#'
#' ## Percentiles interpolados
#' Los tiempos t0, t25, t50 y t75 se calculan mediante interpolación lineal
#' sobre la curva acumulada de germinación. Se aplica un freno biológico
#' adicional: ningún percentil puede ser menor que t0 (primera germinación
#' observada).
#'
#' ## Índice RLG (Relative Light Germination)
#' El RLG se calcula como un valor poblacional único por temperatura,
#' promediando primero las réplicas de luz y oscuridad:
#' \deqn{RLG = \bar{PG}_{luz} / (\bar{PG}_{luz} + \bar{PG}_{oscuridad})}
#'
#' \strong{IMPORTANTE:} RLG NO tiene réplicas estadísticas. Solo debe
#' usarse como métrica descriptiva. NO debe emplearse como variable
#' respuesta en modelos estadísticos (ANOVA, GLM, etc.). Para analizar
#' estadísticamente el efecto de la luz sobre la germinación, se recomienda
#' usar el Porcentaje de Germinación (PG) con un modelo que incluya la
#' interacción Temperatura*Luz.
#'
#' @param data Data frame en formato ancho con columnas de conteos diarios
#'   nombradas como D1, D2, D3... (ej: D0, D1, D2, D3, D4, D5).
#' @param factores Vector de caracteres con los nombres de los factores
#'   explicativos del diseño experimental (ej: \code{c("Temperatura", "Luz")}).
#' @param col_semillas Carácter. Nombre de la columna con el número de
#'   semillas sembradas por réplica. Por defecto \code{"Semillas"}.
#' @param factor_luz Carácter o \code{NULL}. Nombre de la columna que representa
#'   el factor de luz (si existe en el diseño experimental). Si es \code{NULL},
#'   no se calculará el RLG. Por defecto es \code{NULL}.
#' @param umbral_min_rlg Numérico. Porcentaje mínimo de germinación
#'   requerido en al menos una condición (luz u oscuridad) para calcular
#'   el RLG. Por defecto \code{50}.
#' @param etiq_luz Carácter. Etiqueta que identifica la condición de luz
#'   en la columna \code{factor_luz}. Por defecto \code{"Luz"}.
#' @param etiq_osc Carácter. Etiqueta que identifica la condición de
#'   oscuridad. Por defecto \code{"Oscuridad"}.
#' @param solo_luz_cinetica Lógico. Si es \code{TRUE}, restringe el cálculo
#'   de percentiles e índices cinéticos exclusivamente a la condición de
#'   Luz. Si es \code{FALSE}, calcula las métricas para ambas condiciones.
#'   Por defecto es \code{TRUE}.
#'
#' @return Un data frame enriquecido con las siguientes nuevas columnas:
#' \describe{
#'   \item{Total_Germinadas, Exitos, Fracasos}{Conteos binomiales}
#'   \item{PG}{Porcentaje de Germinación}
#'   \item{t0, t25, t50, t75}{Tiempos percentiles interpolados (días)}
#'   \item{TMG}{Tiempo Medio de Germinación}
#'   \item{IVG}{Índice de Velocidad de Germinación}
#'   \item{CVG}{Coeficiente de Velocidad de Germinación}
#'   \item{VGE}{Velocidad de Germinación (1/TMG)}
#'   \item{UNC}{Incertidumbre (Entropía de Shannon)}
#'   \item{SYN}{Índice de Sincronización}
#'   \item{VAR_TG}{Varianza del Tiempo de Germinación}
#'   \item{Z_INDEX}{Índice Z de Czabator}
#'   \item{GRI}{Índice de Tasa de Germinación}
#'   \item{RLG}{Relative Light Germination (un valor por temperatura, NA en el resto)}
#' }
#'
#' @importFrom cli cli_inform cli_warn cli_abort
#' @export
#'
#' @examples
#' \dontrun{
#' # Ejemplo 1: Experimento con Luz y Temperatura
#' matriz_maestra <- preparar_datos_germinacion(
#'   data = datos,
#'   factores = c("Temperatura", "Luz"),
#'   col_semillas = "Semillas",
#'   factor_luz = "Luz",
#'   umbral_min_rlg = 50,
#'   solo_luz_cinetica = TRUE
#' )
#'
#' # Ejemplo 2: Experimento sin luz (solo Humedad y Sustrato)
#' matriz_maestra <- preparar_datos_germinacion(
#'   data = datos,
#'   factores = c("Humedad", "Sustrato"),
#'   col_semillas = "Semillas",
#'   factor_luz = NULL  # No se calculará RLG
#' )
#' }
preparar_datos_germinacion <- function(
    data,
    factores,
    col_semillas = "Semillas",
    factor_luz = NULL,
    umbral_min_rlg = 50,
    etiq_luz = "Luz",
    etiq_osc = "Oscuridad",
    solo_luz_cinetica = TRUE
) {

  # Guardar los nombres de columna originales para comparar al final
  columnas_originales <- names(data)

  # ----------------------------------------------------------
  # 0. VALIDACIONES BÁSICAS DE ENTRADA
  # ----------------------------------------------------------
  if (!is.data.frame(data)) {
    cli::cli_abort("El argumento {.arg data} debe ser un data frame.")
  }
  if (!col_semillas %in% names(data)) {
    cli::cli_abort("No se encontró la columna de semillas: {.val {col_semillas}}.")
  }

  # ----------------------------------------------------------
  # 1. CONTROL DE CASOS LÍMITE: Semillas <= 0 o No Válidas
  # ----------------------------------------------------------
  suppressWarnings(data[[col_semillas]] <- as.numeric(data[[col_semillas]]))
  data[[col_semillas]][is.na(data[[col_semillas]]) | data[[col_semillas]] <= 0] <- NA_real_

  # ----------------------------------------------------------
  # 2. Identificación y ordenamiento cronológico de tiempos
  # ----------------------------------------------------------
  columnas_dias <- names(data)[grepl("^D[0-9]+", names(data))]

  if (length(columnas_dias) == 0) {
    cli::cli_abort("No se encontraron columnas de tiempo con el formato 'D' seguido de números (ej: D1, D2).")
  }

  dias_numericos <- as.numeric(sub("^D", "", columnas_dias))
  orden_cronologico <- order(dias_numericos)
  columnas_dias <- columnas_dias[orden_cronologico]
  dias_numericos <- dias_numericos[orden_cronologico]

  for (col in columnas_dias) {
    suppressWarnings(data[[col]] <- as.numeric(as.character(data[[col]])))
    data[[col]][is.na(data[[col]])] <- 0
  }

  n_filas <- nrow(data)

  # ----------------------------------------------------------
  # 3. [ALERTA DEFENSIVA]: Detección de datos acumulados sospechosos
  # ----------------------------------------------------------
  muestreo_acumulado <- sapply(seq_len(n_filas), function(i) {
    v <- as.numeric(data[i, columnas_dias])
    if (length(v) < 3) return(FALSE)
    all(diff(v) >= 0, na.rm = TRUE) && max(v, na.rm = TRUE) > 0
  })

  if (mean(muestreo_acumulado, na.rm = TRUE) > 0.7) {
    cli::cli_warn(c(
      "!" = "Posible anomalía en los datos: Se detectó un patrón estrictamente creciente en >70% de las filas.",
      "i" = "¿Tus datos están en formato {.strong ACUMULADO}? El IVG e interpolaciones requieren conteos {.strong DIARIOS netos}."
    ))
  }

  # ----------------------------------------------------------
  # 4. Paso Binomial Inicial y Porcentaje de Germinación (PG)
  # ----------------------------------------------------------
  data$Total_Germinadas <- rowSums(data[columnas_dias], na.rm = TRUE)
  data$Exitos <- data$Total_Germinadas
  data$Fracasos <- data[[col_semillas]] - data$Total_Germinadas
  data$PG <- ifelse(!is.na(data[[col_semillas]]),
                    (data$Total_Germinadas / data[[col_semillas]]) * 100,
                    NA_real_)

  # ----------------------------------------------------------
  # 5. VALIDACIÓN Y CONFIGURACIÓN DE FACTOR_LUZ
  # ----------------------------------------------------------
  calcular_rlg <- FALSE
  luz_interna <- NULL
  etiq_luz_norm <- NULL
  etiq_osc_norm <- NULL

  if (!is.null(factor_luz)) {
    if (factor_luz %in% names(data)) {
      luz_interna <- tolower(trimws(as.character(data[[factor_luz]])))
      etiq_luz_norm <- tolower(trimws(etiq_luz))
      etiq_osc_norm <- tolower(trimws(etiq_osc))

      # Verificar que las etiquetas existan en los datos
      niveles_luz <- unique(luz_interna)
      if (etiq_luz_norm %in% niveles_luz && etiq_osc_norm %in% niveles_luz) {
        calcular_rlg <- TRUE
      } else {
        cli::cli_warn(c(
          "!" = "La columna {.val {factor_luz}} no contiene las etiquetas {.val {etiq_luz}} y {.val {etiq_osc}}.",
          "i" = "No se calculará el RLG."
        ))
      }
    } else {
      cli::cli_warn("La columna {.val {factor_luz}} no existe en los datos. No se calculará el RLG.")
    }
  }

  # ----------------------------------------------------------
  # 6. Bucle Único: Índices + Percentiles con Filtros Estrictos
  # ----------------------------------------------------------

  tmg <- numeric(n_filas); ivg <- numeric(n_filas); cvg <- numeric(n_filas)
  unc <- numeric(n_filas); syn <- numeric(n_filas); var_tg <- numeric(n_filas)
  vge <- numeric(n_filas); z_ind <- numeric(n_filas); gri <- numeric(n_filas)
  t0_rep <- numeric(n_filas); t25_rep <- numeric(n_filas)
  t50_rep <- numeric(n_filas); t75_rep <- numeric(n_filas)

  dias_seguros <- ifelse(dias_numericos == 0, 0.1, dias_numericos)

  filas_omitidas_cinetica <- 0

  for (i in seq_len(n_filas)) {
    n_sem <- as.numeric(data[i, col_semillas])
    total_g <- data$Total_Germinadas[i]

    # Determinar si aplicar freno cinético
    aplicar_freno_cinetico <- FALSE
    if (solo_luz_cinetica && !is.null(luz_interna)) {
      es_oscuridad <- (luz_interna[i] != etiq_luz_norm)
      aplicar_freno_cinetico <- es_oscuridad
    }

    if (aplicar_freno_cinetico || is.na(n_sem) || total_g == 0) {
      tmg[i] <- NA_real_; ivg[i] <- NA_real_; cvg[i] <- NA_real_; unc[i] <- NA_real_
      syn[i] <- NA_real_; var_tg[i] <- NA_real_; vge[i] <- NA_real_; z_ind[i] <- NA_real_
      gri[i] <- NA_real_
      t0_rep[i]  <- NA_real_; t25_rep[i] <- NA_real_
      t50_rep[i] <- NA_real_; t75_rep[i] <- NA_real_
      filas_omitidas_cinetica <- filas_omitidas_cinetica + 1
      next
    }

    conteos <- as.numeric(data[i, columnas_dias])

    t0_rep[i]  <- .interpolar_curva(conteos, dias_numericos, 0, n_sem)
    t25_rep[i] <- .interpolar_curva(conteos, dias_numericos, 25, n_sem)
    t50_rep[i] <- .interpolar_curva(conteos, dias_numericos, 50, n_sem)
    t75_rep[i] <- .interpolar_curva(conteos, dias_numericos, 75, n_sem)

    tmg_calc <- sum(conteos * dias_numericos) / total_g
    tmg[i] <- tmg_calc
    ivg[i] <- sum(conteos / dias_seguros)

    cvg[i] <- if (tmg_calc > 0) (total_g / sum(conteos * dias_numericos)) * 100 else NA_real_
    vge[i] <- if (tmg_calc > 0) 1 / tmg_calc else NA_real_

    f_i <- conteos / total_g
    f_i_log <- ifelse(f_i > 0, log2(f_i), 0)
    unc[i] <- -sum(f_i * f_i_log, na.rm = TRUE)

    comb_n_i <- sum(conteos * (conteos - 1)) / 2
    comb_total <- (total_g * (total_g - 1)) / 2
    syn[i] <- if (comb_total > 0) comb_n_i / comb_total else 0

    var_tg[i] <- if (total_g > 1) {
      sum(conteos * (dias_numericos - tmg_calc)^2) / (total_g - 1)
    } else {
      0
    }

    z_ind[i] <- if (total_g > 1) {
      sum(conteos * (total_g - conteos)) / (total_g * (total_g - 1))
    } else {
      0
    }

    pct_diario <- (conteos / n_sem) * 100
    gri[i] <- sum(pct_diario / dias_seguros, na.rm = TRUE)
  }

  data$TMG <- tmg; data$IVG <- ivg; data$CVG <- cvg; data$UNC <- unc
  data$SYN <- syn; data$VAR_TG <- var_tg; data$VGE <- vge; data$Z_INDEX <- z_ind
  data$GRI <- gri
  data$t0  <- t0_rep; data$t25 <- t25_rep; data$t50 <- t50_rep; data$t75 <- t75_rep

  data$t25 <- ifelse(!is.na(data$t25) & !is.na(data$t0) & data$t25 < data$t0, data$t0, data$t25)
  data$t50 <- ifelse(!is.na(data$t50) & !is.na(data$t0) & data$t50 < data$t0, data$t0, data$t50)
  data$t75 <- ifelse(!is.na(data$t75) & !is.na(data$t0) & data$t75 < data$t0, data$t0, data$t75)

  # ----------------------------------------------------------
  # 7. RLG POBLACIONAL: UN SOLO VALOR POR TEMPERATURA (si aplica)
  # ----------------------------------------------------------
  data$RLG <- NA_real_

  if (calcular_rlg) {
    # Identificar la columna de temperatura (primera columna de factores que no es luz)
    col_temperatura <- NULL
    for (f in factores) {
      if (f %in% names(data) && f != factor_luz) {
        col_temperatura <- f
        break
      }
    }

    if (!is.null(col_temperatura)) {
      temp_interna <- tolower(trimws(as.character(data[[col_temperatura]])))
      temperaturas <- unique(temp_interna)

      for (t_actual in temperaturas) {
        filas_temperatura <- which(temp_interna == t_actual)

        pg_luz <- data$PG[temp_interna == t_actual & luz_interna == etiq_luz_norm]
        pg_osc <- data$PG[temp_interna == t_actual & luz_interna == etiq_osc_norm]

        pg_luz <- pg_luz[!is.na(pg_luz)]
        pg_osc <- pg_osc[!is.na(pg_osc)]

        if (length(pg_luz) > 0 && length(pg_osc) > 0) {
          mean_luz <- mean(pg_luz)
          mean_osc <- mean(pg_osc)

          if (mean_luz >= umbral_min_rlg || mean_osc >= umbral_min_rlg) {
            rlg_poblacional <- if (mean_luz == 0 && mean_osc == 0) {
              0.5
            } else {
              mean_luz / (mean_luz + mean_osc)
            }

            data$RLG[filas_temperatura[1]] <- rlg_poblacional
          }
        }
      }
    }
  }

  # ----------------------------------------------------------
  # 8. Identificación de Nuevas Variables y Reporte
  # ----------------------------------------------------------
  columnas_finales <- names(data)
  nuevas_columnas <- setdiff(columnas_finales, columnas_originales)

  if (length(nuevas_columnas) > 0) {
    cli::cli_inform(c(
      "v" = "Matriz maestra enriquecida exitosamente.",
      "i" = "{.strong {length(nuevas_columnas)}} nuevas variables generadas.",
      "*" = "Variables: {.field {paste(nuevas_columnas, collapse = ', ')}}"
    ))
  }

  if (filas_omitidas_cinetica > 0) {
    cli::cli_inform(c(
      "!" = "{.strong {filas_omitidas_cinetica}} de {n_filas} filas fueron omitidas en índices cinéticos.",
      "i" = "Razón: {.val {if (solo_luz_cinetica) 'solo_luz_cinetica = TRUE' else 'semillas = 0 o sin germinación'}}."
    ))
  }

  if (!calcular_rlg) {
    cli::cli_inform(c(
      "i" = "RLG no calculado: {.val {if (is.null(factor_luz)) 'factor_luz = NULL' else 'etiquetas no encontradas'}}."
    ))
  }

  return(data)
}

# ==============================================================================
# FUNCIÓN AUXILIAR INTERNA (NO EXPORTADA)
# ==============================================================================
.interpolar_curva <- function(conteo_num, dias, percentil, n_total) {
  if (is.na(n_total) || n_total <= 0) return(NA_real_)

  acumulado <- cumsum(conteo_num)
  max_germinacion <- max(acumulado, na.rm = TRUE)

  if (is.na(max_germinacion) || max_germinacion == 0) return(NA_real_)

  if (percentil == 0) {
    idx_primer_exito <- which(acumulado > 0)[1]
    return(ifelse(is.na(idx_primer_exito), NA_real_, dias[idx_primer_exito]))
  }

  pct_real_alcanzado <- (max_germinacion / n_total) * 100
  if (pct_real_alcanzado < percentil) return(NA_real_)

  objetivo_semillas <- n_total * (percentil / 100)
  idx_superior <- which(acumulado >= objetivo_semillas)[1]

  if (is.na(idx_superior)) return(NA_real_)
  if (idx_superior == 1) return(dias[1])

  idx_inferior <- idx_superior - 1
  t_inf <- dias[idx_inferior]
  t_sup <- dias[idx_superior]
  n_inf <- acumulado[idx_inferior]
  n_sup <- acumulado[idx_superior]

  if (n_sup == n_inf) return(t_inf)

  return(t_inf + ((objetivo_semillas - n_inf) * (t_sup - t_inf)) / (n_sup - n_inf))
}
