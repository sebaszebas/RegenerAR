#' Ajustar Modelo Estadístico Unificado y Blindado para Germinación
#'
#' @description
#' Ajusta de forma inteligente la estructura del modelo estadístico (Un factor,
#' Interacción o Aditivo) evaluando automáticamente los supuestos biológicos
#' y estadísticos. Se autoconfigura aplicando un GLM Binomial para germinación final,
#' un GLM Gamma para tiempos/percentiles, y Modelos Lineales clásicos con
#' transformación Arcoseno para proporciones y sincronía.
#'
#' @details
#' ## Convención de variable respuesta
#' Cuando \code{variable = "germinacion"}, la función usa internamente
#' \code{cbind(Exitos, Fracasos)} como variable respuesta binomial. Para cualquier
#' otro valor (ej: "TMG", "t50", "SYN"), se usa la columna especificada directamente.
#'
#' ## Estrategias automáticas
#' \itemize{
#'   \item \code{"germinacion"} → GLM Binomial (o Firth si \code{firth = TRUE})
#'   \item \code{TMG, IVG, t0, t25, t50, t75, VAR_TG} → GLM Gamma (link log)
#'   \item \code{SYN, Z_INDEX, RLG} → LM con transformación arcoseno
#'   \item Otras variables → LM clásico
#' }
#'
#' ## Adaptación a diseños asimétricos
#' Si un nivel de "Luz" no tiene datos (común con \code{solo_luz_cinetica = TRUE}),
#' el modelo se adapta automáticamente evaluando solo el factor restante.
#'
#' ## Diagnósticos automáticos
#' Para GLM Binomial, se evalúa sobredispersión (cambiando a Quasibinomial si es
#' necesario) y separación perfecta (sugiriendo Firth si se detecta).
#'
#' ## Paquetes opcionales
#' \itemize{
#'   \item \code{logistf}: Necesario si \code{firth = TRUE}. Instálalo con
#'     \code{instalar_dependencias_opcionales()} o \code{install.packages("logistf")}.
#'   \item \code{performance}: Necesario para el test de sobredispersión.
#'     Instálalo con \code{instalar_dependencias_opcionales()}.
#' }
#'
#' @param data Data frame unificado generado por \code{preparar_datos_germinacion()}.
#' @param factores Vector de caracteres con 1 o 2 factores explicativos.
#' @param variable Carácter. Variable respuesta. Usar \code{"germinacion"} para
#'   análisis binomial, o el nombre de cualquier columna numérica. Por defecto \code{"germinacion"}.
#' @param firth Lógico. Si \code{TRUE} y \code{variable = "germinacion"}, aplica
#'   regresión logística penalizada de Firth. Por defecto \code{FALSE}.
#' @param alpha Numérico. Nivel de significancia. Por defecto \code{0.05}.
#'
#' @return Un objeto S3 de clase \code{"modelo_germinacion_unificado"} con:
#' \describe{
#'   \item{data}{Data frame limpio usado en el modelo}
#'   \item{modelo}{El modelo ajustado (glm, lm, o logistf)}
#'   \item{formula}{La fórmula FINAL utilizada (no la completa inicial)}
#'   \item{factores}{Factores originales solicitados}
#'   \item{factores_activos}{Factores que quedaron tras validaciones}
#'   \item{variable}{Nombre de la variable respuesta}
#'   \item{tipo_modelo}{Tipo: "Un_Factor", "Aditivo", o "Interaccion"}
#'   \item{familia}{Familia final: "GLM_Binomial", "GLM_Quasibinomial", "GLM_Gamma", etc.}
#'   \item{transformacion}{Transformación aplicada a la variable respuesta}
#'   \item{p_interaccion}{P-valor de la interacción (NA si no aplica)}
#'   \item{p_factores}{Vector nombrado con p-valores de cada factor}
#'   \item{alpha}{Nivel de significancia usado}
#'   \item{es_binomial}{Lógico: ¿es análisis binomial?}
#'   \item{es_gamma}{Lógico: ¿es GLM Gamma?}
#'   \item{diagnostico}{Lista con resultados de tests de supuestos}
#' }
#'
#' @importFrom stats glm lm as.formula anova binomial Gamma model.matrix aggregate quasibinomial drop1
#' @importFrom cli cli_abort cli_inform cli_warn
#' @export
#'
#' @examples
#' \dontrun{
#' modelo_tmg <- ajustar_modelo(
#'   data = matriz_maestra,
#'   factores = c("Temperatura", "Luz"),
#'   variable = "TMG"
#' )
#' }
ajustar_modelo <- function(
    data,
    factores,
    variable = "germinacion",
    firth = FALSE,
    alpha = 0.05
) {

  # ============================================================================
  # 1. VALIDACIONES INICIALES
  # ============================================================================
  if (length(factores) < 1 || length(factores) > 2) {
    cli::cli_abort("Debe especificarse 1 o 2 factores explicativos, no {length(factores)}.")
  }

  es_binomial <- (variable == "germinacion")
  columnas_verificar <- if (es_binomial) {
    c("Exitos", "Fracasos", factores)
  } else {
    c(variable, factores)
  }

  faltantes <- setdiff(columnas_verificar, names(data))
  if (length(faltantes) > 0) {
    cli::cli_abort("Faltan las siguientes columnas en el data frame: {.val {faltantes}}")
  }

  data_modelo <- data
  factores_activos <- factores
  diseño_asimetrico_luz <- FALSE

  # ============================================================================
  # 2. ADAPTACIÓN A DISEÑOS ASIMÉTRICOS
  # ============================================================================
  if (!es_binomial) {
    col_luz_candidata <- factores[tolower(trimws(factores)) == "luz"]

    if (length(col_luz_candidata) > 0) {
      niveles_luz <- unique(data_modelo[[col_luz_candidata]])

      for (niv in niveles_luz) {
        valores_nivel <- data_modelo[[variable]][data_modelo[[col_luz_candidata]] == niv]

        if (all(is.na(valores_nivel))) {
          diseño_asimetrico_luz <- TRUE
          factores_activos <- setdiff(factores_activos, col_luz_candidata)

          nivel_con_datos <- niveles_luz[niveles_luz != niv]
          n_filas_antes <- nrow(data_modelo)
          data_modelo <- data_modelo[data_modelo[[col_luz_candidata]] == nivel_con_datos, ]
          n_filas_despues <- nrow(data_modelo)

          cli::cli_inform(c(
            "!" = "Estructura Cinética: Condición {.val {niv}} sin datos para {.val {variable}}.",
            "i" = "Motor adaptado para evaluar {.val {factores_activos}} en condición {.val {nivel_con_datos}}.",
            "v" = "Se eliminaron {n_filas_antes - n_filas_despues} filas sin datos."
          ))
          break
        }
      }
    }

    n_filas_antes <- nrow(data_modelo)
    data_modelo <- data_modelo[!is.na(data_modelo[[variable]]), ]
    n_filas_despues <- nrow(data_modelo)

    if (n_filas_antes != n_filas_despues) {
      cli::cli_inform(c(
        "i" = "Se eliminaron {n_filas_antes - n_filas_despues} filas con NA en {.val {variable}}.",
        "v" = "Datos restantes: {n_filas_despues} observaciones."
      ))
    }

    if (nrow(data_modelo) == 0) {
      cli::cli_abort("No existen datos suficientes para analizar la variable {.val {variable}}.")
    }
  }

  # ============================================================================
  # 3. PREPARACIÓN DE FACTORES
  # ============================================================================
  for (f in factores_activos) {
    if (!is.factor(data_modelo[[f]])) {
      data_modelo[[f]] <- factor(data_modelo[[f]])
    }
  }

  if (!diseño_asimetrico_luz) {
    factores_a_eliminar <- c()
    for (f in factores_activos) {
      if (length(unique(data_modelo[[f]])) < 2) {
        factores_a_eliminar <- c(factores_a_eliminar, f)
      }
    }

    if (length(factores_a_eliminar) > 0) {
      cli::cli_warn(c(
        "!" = "Factores con un único nivel reproducible:",
        "*" = "{.val {factores_a_eliminar}}",
        "i" = "Serán removidos del modelo."
      ))
      factores_activos <- setdiff(factores_activos, factores_a_eliminar)
    }
  }

  if (length(factores_activos) == 0) {
    cli::cli_abort("No existen suficientes niveles reproducibles para ajustar un modelo.")
  }

  # ============================================================================
  # 4. FRENO BIOLÓGICO: Validación de replicación (Ahora con filtrado automático)
  # ============================================================================
  if (!es_binomial) {
    formula_conteo <- stats::as.formula(paste(
      variable, "~", paste(factores_activos, collapse = " + ")
    ))

    conteos_combinacion <- stats::aggregate(
      formula_conteo, data = data_modelo, FUN = length
    )
    names(conteos_combinacion)[ncol(conteos_combinacion)] <- "n_replicas"

    min_replicas_requeridas <- 2
    combinaciones_sin_replicas <- sum(conteos_combinacion$n_replicas < min_replicas_requeridas)

    if (combinaciones_sin_replicas > 0) {
      # Identificar qué combinaciones se van a excluir
      filas_excluir <- conteos_combinacion[conteos_combinacion$n_replicas < min_replicas_requeridas, factores_activos, drop = FALSE]
      filas_excluir_texto <- apply(filas_excluir, 1, function(f) paste(names(f), "=", f, collapse = " & "))

      cli::cli_warn(c(
        "!" = "[Freno Biológico RegenerAR]: Pérdida crítica de datos en algunas combinaciones.",
        "i" = "Variable: {.val {variable}}",
        "i" = "Combinaciones con < {min_replicas_requeridas} réplicas: {length(filas_excluir_texto)}",
        "*" = "{.val {filas_excluir_texto}}",
        "i" = "Estas combinaciones serán excluidas del modelo estadístico."
      ))

      # Filtrar los datos para excluir esas combinaciones
      # Creamos una clave única para filtrar
      data_modelo$clave_filtro <- apply(data_modelo[factores_activos], 1, paste, collapse = "_")
      filas_excluir$clave_filtro <- apply(filas_excluir, 1, paste, collapse = "_")

      data_modelo <- data_modelo[!data_modelo$clave_filtro %in% filas_excluir$clave_filtro, ]
      data_modelo$clave_filtro <- NULL # Limpiar columna temporal

      # Eliminar niveles vacíos de los factores
      for (f in factores_activos) {
        data_modelo[[f]] <- droplevels(data_modelo[[f]])
      }

      cli::cli_inform(c(
        "v" = "Datos filtrados. Observaciones restantes para el modelo: {nrow(data_modelo)}."
      ))
    }

    if (nrow(data_modelo) < (length(factores_activos) + 2)) {
      cli::cli_abort(c(
        "x" = "[Freno Biológico RegenerAR]: Pérdida crítica de datos total.",
        "i" = "Tras excluir combinaciones problemáticas, no quedan suficientes datos para ajustar el modelo.",
        "i" = "Observaciones totales restantes: {nrow(data_modelo)}"
      ))
    }
  }

  # ============================================================================
  # 5. DETERMINACIÓN DE ESTRATEGIA
  # ============================================================================
  estrategia <- "LM_Clasico"
  transformacion <- "Ninguna"
  diagnostico <- list()

  if (es_binomial) {
    estrategia <- if (firth) "GLM_Firth" else "GLM_Binomial"
  } else if (variable %in% c("TMG", "IVG", "t0", "t25", "t50", "t75", "VAR_TG")) {
    estrategia <- "GLM_Gamma"
    if (any(data_modelo[[variable]] <= 0)) {
      data_modelo[[variable]] <- ifelse(
        data_modelo[[variable]] <= 0,
        data_modelo[[variable]] + 0.01,
        data_modelo[[variable]]
      )
      transformacion <- "Ajuste de coto inferior (+0.01) para convergencia Gamma"
      cli::cli_warn(c(
        "!" = "Se detectaron valores ≤ 0 en {.val {variable}}.",
        "i" = "Se aplicó ajuste de +0.01 para permitir convergencia del GLM Gamma.",
        "i" = "Interprete los resultados con precaución."
      ))
    }
  } else if (variable %in% c("SYN", "Z_INDEX", "RLG")) {
    valores_var <- data_modelo[[variable]]
    if (any(valores_var < 0 | valores_var > 1, na.rm = TRUE)) {
      cli::cli_abort(c(
        "x" = "La variable {.val {variable}} contiene valores fuera del rango [0, 1].",
        "i" = "La transformación arcoseno requiere proporciones en el intervalo [0, 1].",
        "i" = "Rango detectado: [{min(valores_var, na.rm = TRUE)}, {max(valores_var, na.rm = TRUE)}]"
      ))
    }
    estrategia <- "LM_Arcoseno"
    data_modelo[[variable]] <- asin(sqrt(data_modelo[[variable]]))
    transformacion <- "Arcoseno de la raíz cuadrada: asin(sqrt(x))"
  }

  # ============================================================================
  # 6. VERIFICACIÓN DE PAQUETES OPCIONALES (logistf, performance)
  # ============================================================================
  if (estrategia == "GLM_Firth" && !requireNamespace("logistf", quietly = TRUE)) {
    cli::cli_abort(c(
      "x" = "El paquete {.pkg logistf} es necesario para la regresión de Firth.",
      "i" = "Instálalo con: {.code install.packages('logistf')}",
      "i" = "O ejecuta: {.code RegenerAR::instalar_dependencias_opcionales()}"
    ))
  }

  # ============================================================================
  # 7. AJUSTE DEL MODELO COMPLETO (con interacción)
  # ============================================================================
  formula_completa <- stats::as.formula(paste(
    if (es_binomial) "cbind(Exitos, Fracasos)" else variable,
    "~",
    paste(factores_activos, collapse = " * ")
  ))

  modelo_completo <- .ajustar_motor(formula_completa, data_modelo, estrategia)

  # Registrar si Firth falló (para ajustar mensajes posteriores)
  firth_fallo <- FALSE

  if (estrategia == "GLM_Firth" && inherits(modelo_completo, "try-error")) {
    firth_fallo <- TRUE
    cli::cli_warn(c(
      "!" = "Firth no convergió numéricamente.",
      "i" = "Conmutando a GLM Binomial Clásico de emergencia."
    ))
    estrategia <- "GLM_Binomial"
    modelo_completo <- .ajustar_motor(formula_completa, data_modelo, estrategia)
  }

  # ============================================================================
  # 8. SELECCIÓN DE ESTRUCTURA FINAL
  # ============================================================================
  p_interaccion <- NA_real_

  if (length(factores_activos) == 1) {
    modelo_final <- modelo_completo
    formula_final <- formula_completa
    tipo_modelo <- "Un_Factor"
  } else {
    formula_aditiva <- stats::as.formula(paste(
      if (es_binomial) "cbind(Exitos, Fracasos)" else variable,
      "~",
      paste(factores_activos, collapse = " + ")
    ))
    modelo_aditivo <- .ajustar_motor(formula_aditiva, data_modelo, estrategia)

    if (estrategia == "GLM_Firth") {
      comparacion <- stats::anova(modelo_aditivo, modelo_completo)
      p_interaccion <- comparacion$P[2]
    } else {
      test_type <- if (estrategia %in% c("GLM_Binomial", "GLM_Gamma")) "Chisq" else NULL
      comparacion <- stats::anova(modelo_aditivo, modelo_completo, test = test_type)
      idx_p <- grep("Pr\\(>", names(comparacion))
      p_interaccion <- if (length(idx_p) > 0) comparacion[2, idx_p] else NA_real_
    }

    if (!is.na(p_interaccion) && p_interaccion < alpha) {
      modelo_final <- modelo_completo
      formula_final <- formula_completa
      tipo_modelo <- "Interaccion"
      cli::cli_inform(c(
        "v" = "Interacción significativa (p = {format.pval(p_interaccion, digits = 3)}).",
        "i" = "Estructura final: {.strong Interacción}"
      ))
    } else {
      modelo_final <- modelo_aditivo
      formula_final <- formula_aditiva
      tipo_modelo <- "Aditivo"
      cli::cli_inform(c(
        "i" = "Interacción no significativa (p = {format.pval(p_interaccion, digits = 3)}).",
        "v" = "Estructura final: {.strong Aditivo}"
      ))
    }
  }

  # ============================================================================
  # 9. DIAGNÓSTICOS DE SUPUESTOS
  # ============================================================================
  familia_reporte <- estrategia

  if (estrategia == "GLM_Binomial") {
    # 🔧 FIX: Envolver TODO el bloque de diagnóstico para capturar warnings de cbind
    # en cualquier parte (detectseparation, aggregate, etc.)

    withCallingHandlers({

      # 9a. Test de sobredispersión (requiere 'performance')
      if (requireNamespace("performance", quietly = TRUE)) {
        diag_obj <- performance::check_overdispersion(modelo_final, verbose = FALSE)
        p_disp <- diag_obj$p_value
        ratio_disp <- diag_obj$dispersion_ratio

        diagnostico$dispersion <- list(
          ratio = ratio_disp,
          p_value = p_disp,
          severa = !is.na(p_disp) && p_disp < alpha && ratio_disp > 1
        )

        if (diagnostico$dispersion$severa) {
          familia_reporte <- "GLM_Quasibinomial"
          modelo_final <- stats::glm(
            formula_final,
            data = data_modelo,
            family = stats::quasibinomial(link = "logit")
          )
          cli::cli_warn(c(
            "!" = "Sobredispersión crítica detectada.",
            "i" = "Ratio: {round(ratio_disp, 2)} | p = {format.pval(p_disp, digits = 3)}",
            "v" = "Distribución final: {.strong QUASIBINOMIAL}"
          ))
        } else {
          cli::cli_inform(c(
            "v" = "Dispersión correcta (Ratio: {round(ratio_disp, 2)}, p = {format.pval(p_disp, digits = 3)}).",
            "i" = "Distribución final: {.strong BINOMIAL}"
          ))
        }
      } else {
        cli::cli_inform(c(
          "i" = "Test de sobredispersión omitido (paquete {.pkg performance} no instalado).",
          "i" = "Para activarlo: {.code install.packages('performance')}"
        ))
      }

      # 9b. Test de separación perfecta
      matriz_y <- as.matrix(data_modelo[, c("Exitos", "Fracasos")])

      test_sep <- tryCatch(
        detectseparation::detect_separation(
          x = stats::model.matrix(formula_final, data = data_modelo),
          y = matriz_y,
          family = stats::binomial()
        ),
        error = function(e) {
          cli::cli_warn("No se pudo ejecutar el test de separación: {e$message}")
          NULL
        }
      )

      if (!is.null(test_sep) && isTRUE(test_sep$outcome)) {
        columnas_agrupar <- c("Exitos", "Fracasos", factores_activos)
        resumen_sep <- stats::aggregate(
          cbind(Exitos, Fracasos) ~ .,
          data = data_modelo[columnas_agrupar],
          FUN = sum
        )
        resumen_sep$tasa <- resumen_sep$Exitos / (resumen_sep$Exitos + resumen_sep$Fracasos)
        filas_problema <- resumen_sep[resumen_sep$tasa == 0 | resumen_sep$tasa == 1, ]

        if (nrow(filas_problema) > 0) {
          combinaciones_problema <- apply(filas_problema, 1, function(fila) {
            elementos <- sapply(factores_activos, function(fa) paste0(fa, "=", fila[fa]))
            paste0("(", paste(elementos, collapse = " × "), ")")
          })

          diagnostico$separacion <- list(
            detectada = TRUE,
            combinaciones = combinaciones_problema
          )

          if (!firth) {
            cli::cli_warn(c(
              "!" = "Separación perfecta detectada en:",
              "*" = "{.val {combinaciones_problema}}",
              "i" = "Sugerencia: re-ejecute con {.code firth = TRUE}",
              "i" = "Nota: Firth requiere el paquete {.pkg logistf}. Instálalo con:",
              " " = "{.code install.packages('logistf')}"
            ))
          } else {
            cli::cli_warn(c(
              "!" = "Separación perfecta detectada en:",
              "*" = "{.val {combinaciones_problema}}",
              "x" = "Firth fue solicitado pero no convergió numéricamente.",
              "i" = "El modelo final usa GLM Binomial clásico (resultados pueden ser inestables).",
              "i" = "Alternativas:",
              " " = "1. Revisar si hay combinaciones con 0% o 100% de germinación",
              " " = "2. Considerar colapsar niveles de factores",
              " " = "3. Usar un modelo bayesiano con priors informativos"
            ))
          }
        }
      }

    }, warning = function(w) {
      # 🔧 Filtrar SOLO el warning específico de cbind
      if (grepl("number of rows of result is not a multiple", w$message)) {
        invokeRestart("muffleWarning")
      }
      # Otros warnings pasan normalmente
    })

  } else if (estrategia == "GLM_Firth") {
    cli::cli_inform(c(
      "v" = "Distribución final: {.strong Regresión Logística Penalizada de Firth}",
      "i" = "Errores estándar corregidos ante separación perfecta."
    ))
  } else if (estrategia == "GLM_Gamma") {
    cli::cli_inform(c(
      "v" = "Distribución final: {.strong GLM Gamma (link log)}"
    ))
  } else if (estrategia == "LM_Arcoseno") {
    cli::cli_inform(c(
      "v" = "Distribución final: {.strong LM con transformación Arcoseno}",
      "i" = "Nota: Los resultados del post-hoc se back-transformarán automáticamente."
    ))
  } else {
    cli::cli_inform(c(
      "v" = "Distribución final: {.strong LM Clásico (Mínimos Cuadrados)}"
    ))
  }

  # ============================================================================
  # 10. EXTRACCIÓN DE P-VALORES
  # ============================================================================
  p_factores <- NULL

  if (tipo_modelo == "Un_Factor") {
    formula_nula <- stats::as.formula(paste(
      if (es_binomial) "cbind(Exitos, Fracasos)" else variable, "~ 1"
    ))
    modelo_nulo <- .ajustar_motor(formula_nula, data_modelo, estrategia)

    test_type <- if (estrategia == "GLM_Gamma") "F" else if (estrategia %in% c("GLM_Binomial")) "Chisq" else NULL
    comparacion_nulo <- stats::anova(modelo_nulo, modelo_final, test = test_type)

    idx_p <- grep("Pr\\(>", names(comparacion_nulo))
    p_global <- if (length(idx_p) > 0) comparacion_nulo[2, idx_p] else NA_real_
    p_factores <- stats::setNames(p_global, factores_activos)

    sig_texto <- if (!is.na(p_global) && p_global < alpha) "SIGNIFICATIVO" else "NO SIGNIFICATIVO"
    cli::cli_inform(c(
      "i" = "Variable [{.val {variable}}] | Factor {.val {factores_activos}}: {.strong {sig_texto}} (p = {format.pval(p_global, digits = 3)})."
    ))

  } else if (tipo_modelo == "Interaccion") {
    if (familia_reporte == "GLM_Firth") {
      anova_int <- stats::anova(modelo_final)
      p_factores <- anova_int$P
      names(p_factores) <- rownames(anova_int)
    } else {
      test_type <- if (familia_reporte %in% c("GLM_Binomial", "GLM_Quasibinomial", "GLM_Gamma")) "Chisq" else "F"
      anova_int <- stats::anova(modelo_final, test = test_type)
      idx_p <- grep("Pr\\(>", names(anova_int))
      if (length(idx_p) > 0) {
        p_factores <- anova_int[[idx_p]]
        names(p_factores) <- rownames(anova_int)
        p_factores <- p_factores[!is.na(p_factores)]
      }
    }

    termino_interaccion <- paste(factores_activos, collapse = ":")
    p_int_real <- p_factores[termino_interaccion]

    if (length(p_int_real) > 0 && !is.na(p_int_real)) {
      sig_texto <- if (p_int_real < alpha) "SIGNIFICATIVO" else "NO SIGNIFICATIVO"
      cli::cli_inform(c(
        "i" = "Variable [{.val {variable}}] | Interacción {.val {termino_interaccion}}: {.strong {sig_texto}} (p = {format.pval(p_int_real, digits = 3)})."
      ))
    }

  } else if (tipo_modelo == "Aditivo") {
    if (familia_reporte == "GLM_Firth") {
      anova_adit <- stats::anova(modelo_final)
      p_factores <- anova_adit$P
      names(p_factores) <- rownames(anova_adit)
    } else {
      test_type <- if (familia_reporte %in% c("GLM_Binomial", "GLM_Gamma")) "Chisq" else "F"
      anova_adit <- stats::drop1(modelo_final, test = test_type)
      col_p <- intersect(c("Pr(>Chi)", "Pr(>F)"), names(anova_adit))

      if (length(col_p) > 0) {
        p_raw <- anova_adit[[col_p]]
        nombres_raw <- rownames(anova_adit)
        filas_validas <- nombres_raw != "<none>"
        p_factores <- p_raw[filas_validas]
        names(p_factores) <- nombres_raw[filas_validas]
      }
    }

    for (f in factores_activos) {
      p <- p_factores[f]
      sig_texto <- if (!is.na(p) && p < alpha) "SIGNIFICATIVO" else "NO SIGNIFICATIVO"
      cli::cli_inform(c(
        "i" = "Variable [{.val {variable}}] | Factor {.val {f}}: {.strong {sig_texto}} (p = {format.pval(p, digits = 3)})."
      ))
    }
  }

  # ============================================================================
  # 11. OBJETO DE SALIDA S3
  # ============================================================================
  resultado <- list(
    data             = data_modelo,
    modelo           = modelo_final,
    formula          = formula_final,
    factores         = factores,
    factores_activos = factores_activos,
    variable         = variable,
    tipo_modelo      = tipo_modelo,
    familia          = familia_reporte,
    transformacion   = transformacion,
    p_interaccion    = p_interaccion,
    p_factores       = p_factores,
    alpha            = alpha,
    es_binomial      = es_binomial,
    es_gamma         = (estrategia == "GLM_Gamma"),
    diagnostico      = diagnostico
  )

  class(resultado) <- "modelo_germinacion_unificado"

  cli::cli_inform(c(
    "v" = "Modelo {.strong {tipo_modelo}} ajustado exitosamente.",
    "i" = "Familia: {.val {familia_reporte}} | Fórmula: {.code {deparse(formula_final)}}"
  ))

  return(resultado)
}

# ==============================================================================
# FUNCIÓN AUXILIAR INTERNA (NO EXPORTADA)
# ==============================================================================
.ajustar_motor <- function(form, datos, est) {
  tryCatch(
    withCallingHandlers({
      switch(est,
             "GLM_Binomial" = stats::glm(form, data = datos, family = stats::binomial(link = "logit")),
             "GLM_Firth"    = {
               if (!requireNamespace("logistf", quietly = TRUE)) {
                 cli::cli_abort(c(
                   "x" = "El paquete {.pkg logistf} es necesario.",
                   "i" = "Instálalo con: {.code install.packages('logistf')}"
                 ))
               }
               logistf::logistf(form, data = datos)
             },
             "GLM_Gamma"    = stats::glm(form, data = datos, family = stats::Gamma(link = "log")),
             "LM_Arcoseno"  = stats::lm(form, data = datos),
             "LM_Clasico"   = stats::lm(form, data = datos)
      )
    }, warning = function(w) {
      # 🔧 Filtrar SOLO el warning específico de cbind en glm binomial
      if (grepl("number of rows of result is not a multiple", w$message)) {
        invokeRestart("muffleWarning")
      }
      # Otros warnings pasan normalmente
    }),
    error = function(e) {
      if (est == "GLM_Firth") {
        structure(list(message = e$message), class = "try-error")
      } else {
        cli::cli_abort("Error al ajustar el modelo {.val {est}}: {e$message}")
      }
    }
  )
}
