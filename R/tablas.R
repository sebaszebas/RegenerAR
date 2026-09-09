#' Generar Tablas Resumen Científicas (Formato de Publicación)
#'
#' @description
#' Construye reportes tabulares con alineación simétrica y estética para la consola de R.
#' Soporta tablas de significancia individuales, matrices cruzadas de percentiles (t0 a t75),
#' reportes del Índice de Fotoblastismo (RLG) y la tabla maestra de óptimos (Germinación + TMG + IVG).
#'
#' @details
#' ## Modos de operación
#' \itemize{
#'   \item Si \code{objeto} es un modelo ajustado → delega a \code{post_hoc_germinacion()}
#'   \item Si \code{objeto} es un data frame + \code{tipo = "percentiles"} → matriz de t0-t75
#'   \item Si \code{tipo = "fotoblastismo"} → reporte de RLG por temperatura
#'   \item Si \code{tipo = "optimo"} → tabla maestra con Germinación + TMG + IVG fusionados
#' }
#'
#' ## Nota sobre el formato
#' Los valores se formatean con centrado estético para visualización en consola.
#' Los valores numéricos quedan como caracteres con espacios. Si necesita exportar
#' a Excel/Word para análisis posterior, use \code{limpiar_tabla()} para remover
#' los espacios de relleno.
#'
#' @param objeto Data frame maestro generado por \code{preparar_datos_germinacion()}
#'   o un objeto de clase \code{"modelo_germinacion_unificado"}.
#' @param tipo Carácter. Tipo de reporte: \code{"percentiles"}, \code{"fotoblastismo"}
#'   o \code{"optimo"}. Se ignora si \code{objeto} es un modelo ajustado.
#' @param factores Vector de caracteres con los nombres de los factores explicativos.
#'   Si es \code{NULL}, se autodetectan (máximo 2). Por defecto \code{NULL}.
#' @param alpha Numérico. Nivel de significancia. Por defecto 0.05.
#' @param firth Lógico. Si \code{TRUE} y \code{tipo = "optimo"}, aplica regresión
#'   de Firth para el submodelo de germinación final. Por defecto \code{FALSE}.
#' @param verbose Lógico. Si \code{FALSE}, silencia los mensajes de los modelos
#'   internos en \code{tipo = "optimo"}. Por defecto \code{TRUE}.
#'
#' @return Un data frame formateado con centrado estético para visualización en consola.
#'
#' @importFrom dplyr group_by across all_of summarise filter left_join select n
#' @importFrom cli cli_abort cli_inform cli_warn
#' @export
#'
#' @examples
#' \dontrun{
#' # Tabla de percentiles
#' tabla_perc <- reportar_resultados(matriz_maestra, tipo = "percentiles",
#'                                    factores = c("Temperatura", "Luz"))
#'
#' # Tabla maestra de óptimos
#' tabla_opt <- reportar_resultados(matriz_maestra, tipo = "optimo",
#'                                   factores = c("Temperatura", "Luz"))
#'
#' # Post-hoc directo desde un modelo
#' posthoc <- reportar_resultados(modelo_ajustado)
#' }
reportar_resultados <- function(objeto, tipo = "percentiles", factores = NULL,
                                alpha = 0.05, firth = FALSE, verbose = TRUE) {

  # ============================================================================
  # VALIDACIONES INICIALES
  # ============================================================================
  # CASO A: El usuario introduce un modelo → delegar a post_hoc
  if (inherits(objeto, "modelo_germinacion_unificado")) {
    return(post_hoc_germinacion(objeto, alpha = alpha))
  }

  if (!is.data.frame(objeto)) {
    cli::cli_abort(c(
      "x" = "El objeto debe ser la matriz maestra o un modelo de RegenerAR.",
      "i" = "Clase detectada: {.cls {class(objeto)}}"
    ))
  }

  if (!tipo %in% c("percentiles", "fotoblastismo", "optimo")) {
    cli::cli_abort(c(
      "x" = "Tipo de reporte no reconocido: {.val {tipo}}",
      "i" = "Opciones válidas: {.val percentiles}, {.val fotoblastismo}, {.val optimo}"
    ))
  }

  # ============================================================================
  # DETECCIÓN AUTOMÁTICA DE FACTORES
  # ============================================================================
  if (is.null(factores)) {
    columnas_dias <- names(objeto)[grepl("^D[0-9]+$", names(objeto))]
    columnas_excluir <- c(
      columnas_dias,
      "Total_Germinadas", "Exitos", "Fracasos", "PG", "TMG", "IVG", "CVG",
      "UNC", "SYN", "VAR_TG", "VGE", "Z_INDEX", "GRI",
      "t0", "t25", "t50", "t75", "RLG",
      "Rep", "Replica", "Replicas", "Placa", "ID", "id", "rep", "replica"
    )

    factores_potenciales <- setdiff(names(objeto), columnas_excluir)
    factores_detectados <- head(factores_potenciales, 2)

    if (length(factores_detectados) == 0) {
      cli::cli_abort(c(
        "x" = "No se pudieron autodetectar los factores explicativos.",
        "i" = "Especifique manualmente: {.code factores = c('Factor1', 'Factor2')}"
      ))
    }

    if (length(factores_potenciales) > 2) {
      cli::cli_warn(c(
        "!" = "Se detectaron {length(factores_potenciales)} columnas candidatas a factores.",
        "i" = "Usando las primeras 2: {.val {factores_detectados}}",
        "i" = "Si no son correctas, especifique: {.code factores = c(...)}"
      ))
    }
  } else {
    factores_detectados <- factores
  }

  if (length(factores_detectados) < 1 || length(factores_detectados) > 2) {
    cli::cli_abort(c(
      "x" = "Se detectaron {length(factores_detectados)} factores: {.val {factores_detectados}}",
      "i" = "RegenerAR solo admite 1 o 2 factores explicativos.",
      "i" = "Especifique manualmente: {.code factores = c('Factor1', 'Factor2')}"
    ))
  }

  # Limpieza de espacios invisibles
  for (f_col in factores_detectados) {
    if (f_col %in% names(objeto)) {
      objeto[[f_col]] <- trimws(as.character(objeto[[f_col]]))
    }
  }

  # ============================================================================
  # CASO B: TABLA MAESTRA DE ÓPTIMOS (Germinación + TMG + IVG)
  # ============================================================================
  if (tipo == "optimo") {
    if (verbose) {
      cli::cli_inform(c(
        "i" = "Generando Tabla Maestra de Óptimos Fisiológicos...",
        "i" = "Ajustando 3 modelos internos: Germinación, TMG, IVG"
      ))
    }

    # Helper para ejecutar modelos silenciosamente
    .ejecutar_modelo_silencioso <- function(expr) {
      if (verbose) expr else suppressMessages(suppressWarnings(expr))
    }

    # 🔧 MOTOR INTELIGENTE: Exclusión, Modelado y Fallback Automático
    .procesar_variable_cinetica <- function(data, factores, variable) {
      # 1. Contar réplicas válidas (no-NA) por combinación de factores
      conteos <- data |>
        dplyr::group_by(dplyr::across(dplyr::all_of(factores))) |>
        dplyr::summarise(n_validos = sum(!is.na(.data[[variable]])), .groups = "drop")

      niveles_excluidos <- conteos[conteos$n_validos < 2, factores, drop = FALSE]
      niveles_validos <- conteos[conteos$n_validos >= 2, factores, drop = FALSE]

      # Si no hay niveles válidos, todo es guión
      if (nrow(niveles_validos) == 0) {
        out <- niveles_excluidos
        out[[variable]] <- "-"
        out[[paste0(variable, "_let")]] <- "-"
        return(out)
      }

      # 2. Filtrar datos para el modelo
      data_filtrada <- dplyr::inner_join(data, niveles_validos, by = factores)

      # 3. Ajustar modelo (o fallback automático a Media ± DE si falla)
      tabla_modelo <- tryCatch({
        mod <- .ejecutar_modelo_silencioso(ajustar_modelo(data_filtrada, factores, variable, alpha = alpha))
        res <- .ejecutar_modelo_silencioso(post_hoc_germinacion(mod, alpha = alpha))
        names(res) <- trimws(names(res))
        col_let <- names(res)[grepl("_let", names(res))]
        if (length(col_let) > 0) names(res)[names(res) == col_let] <- paste0(variable, "_let")
        res
      }, error = function(e) {
        if (verbose) cli::cli_warn("Modelo {variable} falló. Mostrando Media ± DE automáticamente.")

        # Fallback descriptivo automático
        data_clean <- data_filtrada[!is.na(data_filtrada[[variable]]), ]
        resumen <- data_clean |>
          dplyr::group_by(dplyr::across(dplyr::all_of(factores))) |>
          dplyr::summarise(media = mean(.data[[variable]], na.rm = TRUE),
                           sd = stats::sd(.data[[variable]], na.rm = TRUE),
                           n = dplyr::n(), .groups = "drop")
        out <- data.frame(resumen[factores], check.names = FALSE)
        out[[variable]] <- vapply(seq_len(nrow(resumen)), function(i) {
          m <- resumen$media[i]; s <- resumen$sd[i]; n <- resumen$n[i]
          if (is.na(s) || s == 0 || n < 2) return(as.character(round(m, 1)))
          return(paste0(round(m, 1), " ± ", round(s, 1)))
        }, character(1))
        out[[paste0(variable, "_let")]] <- "-"
        out
      })

      # 4. Reincorporar niveles excluidos con guión "-"
      if (nrow(niveles_excluidos) > 0) {
        tabla_excluidos <- niveles_excluidos
        tabla_excluidos[[variable]] <- "-"
        tabla_excluidos[[paste0(variable, "_let")]] <- "-"
        tabla_modelo <- dplyr::bind_rows(tabla_modelo, tabla_excluidos)
      }
      return(tabla_modelo)
    }

    # 1. Post-hoc de Germinación Final
    mod_pg <- .ejecutar_modelo_silencioso(
      ajustar_modelo(objeto, factores = factores_detectados,
                     variable = "germinacion", alpha = alpha, firth = firth)
    )
    tabla_pg_letras <- .ejecutar_modelo_silencioso(
      post_hoc_germinacion(mod_pg, alpha = alpha)
    )
    names(tabla_pg_letras) <- trimws(names(tabla_pg_letras))

    # Control de luz para variables cinéticas
    factores_cinetica <- factores_detectados
    col_luz_candidata <- factores_detectados[tolower(trimws(factores_detectados)) == "luz"]

    if (length(col_luz_candidata) > 0 && "TMG" %in% names(objeto)) {
      niveles_luz <- unique(objeto[[col_luz_candidata]])
      for (niv in niveles_luz) {
        if (all(is.na(objeto$TMG[objeto[[col_luz_candidata]] == niv]))) {
          factores_cinetica <- setdiff(factores_cinetica, col_luz_candidata)
          break
        }
      }
    }

    # 2. Procesar TMG y IVG con el motor inteligente
    tabla_tmg_letras <- .procesar_variable_cinetica(objeto, factores_cinetica, "TMG")
    tabla_ivg_letras <- .procesar_variable_cinetica(objeto, factores_cinetica, "IVG")

    # 3. FUSIÓN DE LAS TRES TABLAS
    tabla_fusionada <- tabla_pg_letras

    # A. Acoplar TMG
    if (!is.null(tabla_tmg_letras)) {
      col_valor_tmg <- names(tabla_tmg_letras)[grepl("TMG", names(tabla_tmg_letras)) & !grepl("_let", names(tabla_tmg_letras))]
      by_tmg <- intersect(names(tabla_fusionada), intersect(names(tabla_tmg_letras), factores_cinetica))

      tabla_fusionada <- dplyr::left_join(
        tabla_fusionada,
        tabla_tmg_letras[, c(by_tmg, col_valor_tmg, "TMG_let"), drop = FALSE],
        by = by_tmg
      )

      if (length(col_luz_candidata) > 0 && col_luz_candidata %in% names(tabla_fusionada)) {
        es_oscuridad <- tolower(trimws(as.character(tabla_fusionada[[col_luz_candidata]]))) != "luz"
        tabla_fusionada[[col_valor_tmg]][es_oscuridad] <- "-"
        tabla_fusionada[["TMG_let"]][es_oscuridad] <- "-"
      }
      titulo_final_tmg <- col_valor_tmg
    } else {
      titulo_final_tmg <- "TMG"
      tabla_fusionada[[titulo_final_tmg]] <- "-"
      tabla_fusionada[["TMG_let"]] <- "-"
    }

    # B. Acoplar IVG
    if (!is.null(tabla_ivg_letras)) {
      col_valor_ivg <- names(tabla_ivg_letras)[grepl("IVG", names(tabla_ivg_letras)) & !grepl("_let", names(tabla_ivg_letras))]
      by_ivg <- intersect(names(tabla_fusionada), intersect(names(tabla_ivg_letras), factores_cinetica))

      tabla_fusionada <- dplyr::left_join(
        tabla_fusionada,
        tabla_ivg_letras[, c(by_ivg, col_valor_ivg, "IVG_let"), drop = FALSE],
        by = by_ivg
      )

      if (length(col_luz_candidata) > 0 && col_luz_candidata %in% names(tabla_fusionada)) {
        es_oscuridad <- tolower(trimws(as.character(tabla_fusionada[[col_luz_candidata]]))) != "luz"
        tabla_fusionada[[col_valor_ivg]][es_oscuridad] <- "-"
        tabla_fusionada[["IVG_let"]][es_oscuridad] <- "-"
      }
      titulo_final_ivg <- col_valor_ivg
    } else {
      titulo_final_ivg <- "IVG"
      tabla_fusionada[[titulo_final_ivg]] <- "-"
      tabla_fusionada[["IVG_let"]] <- "-"
    }

    # 4. Centrado estético final
    col_germinacion <- names(tabla_fusionada)[grepl("Germinaci", names(tabla_fusionada))]
    col_let_pg <- names(tabla_fusionada)[grepl("_let", names(tabla_fusionada)) &
                                           !names(tabla_fusionada) %in% c("TMG_let", "IVG_let")]

    cols_a_formatear <- c(col_germinacion, col_let_pg, titulo_final_tmg, "TMG_let",
                          titulo_final_ivg, "IVG_let")

    for (col_f in cols_a_formatear) {
      if (col_f %in% names(tabla_fusionada)) {
        tabla_fusionada <- .centrar_columna_texto(tabla_fusionada, col_f, col_f)
      }
    }

    cli::cli_inform(c(
      "v" = "Tabla maestra generada exitosamente.",
      "i" = "Óptimo fisiológico: maximizar Germinación e IVG, minimizar TMG."
    ))

    return(as.data.frame(tabla_fusionada))
  }

  # ============================================================================
  # CASO C: REPORTE DE PERCENTILES (t0, t25, t50, t75)
  # ============================================================================
  if (tipo == "percentiles") {
    if (!all(c("t0", "t25", "t50", "t75") %in% names(objeto))) {
      cli::cli_abort("La matriz no contiene las columnas de percentiles (t0, t25, t50, t75).")
    }

    resumen <- objeto |>
      dplyr::group_by(dplyr::across(dplyr::all_of(factores_detectados))) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(c("t0", "t25", "t50", "t75")),
          list(
            m = function(x) mean(x, na.rm = TRUE),
            s = function(x) stats::sd(x, na.rm = TRUE)
          ),
          .names = "{col}_{fn}"
        ),
        n_rep = dplyr::n(),
        .groups = "drop"
      )

    # 🔧 BLINDAJE DE MONOTONÍA BIOLÓGICA
    # Si por falta de réplicas (n=1) la media de un percentil superior es menor
    # que el inferior, forzamos la igualdad para mantener coherencia fisiológica.
    resumen$t75_m <- ifelse(resumen$t75_m < resumen$t50_m, resumen$t50_m, resumen$t75_m)
    resumen$t50_m <- ifelse(resumen$t50_m < resumen$t25_m, resumen$t25_m, resumen$t50_m)
    resumen$t25_m <- ifelse(resumen$t25_m < resumen$t0_m,  resumen$t0_m,  resumen$t25_m)

    fmt_celda <- function(m, s, n) {
      if (is.na(m) || is.nan(m)) return("-")
      # Si hay 1 sola réplica, el sd es NA. Mostramos solo el valor.
      if (is.na(s) || s == 0 || n < 2) return(as.character(round(m, 1)))
      return(paste0(round(m, 1), " ± ", round(s, 1)))
    }

    tabla_fmt <- data.frame(resumen[factores_detectados], check.names = FALSE)

    # Usamos los valores YA corregidos por el blindaje biológico
    tabla_fmt$t0  <- vapply(seq_len(nrow(resumen)), function(i) fmt_celda(resumen$t0_m[i],  resumen$t0_s[i],  resumen$n_rep[i]), character(1))
    tabla_fmt$t25 <- vapply(seq_len(nrow(resumen)), function(i) fmt_celda(resumen$t25_m[i], resumen$t25_s[i], resumen$n_rep[i]), character(1))
    tabla_fmt$t50 <- vapply(seq_len(nrow(resumen)), function(i) fmt_celda(resumen$t50_m[i], resumen$t50_s[i], resumen$n_rep[i]), character(1))
    tabla_fmt$t75 <- vapply(seq_len(nrow(resumen)), function(i) fmt_celda(resumen$t75_m[i], resumen$t75_s[i], resumen$n_rep[i]), character(1))

    for (perc in c("t0", "t25", "t50", "t75")) {
      tabla_fmt <- .centrar_columna_texto(tabla_fmt, perc, perc)
    }

    cli::cli_inform(c(
      "v" = "Tabla de percentiles generada.",
      "i" = "Valores: Media ± DE (días). Se aplicó blindaje de monotonía biológica."
    ))

    return(as.data.frame(tabla_fmt))
  }
  # ============================================================================
  # CASO D: REPORTE DE FOTOBLASTISMO (RLG)
  # ============================================================================
  if (tipo == "fotoblastismo") {
    if (!"RLG" %in% names(objeto)) {
      cli::cli_abort("La matriz no contiene la columna de Fotoblastismo (RLG).")
    }

    factores_sin_luz <- factores_detectados[tolower(trimws(factores_detectados)) != "luz"]
    if (length(factores_sin_luz) == 0) factores_sin_luz <- factores_detectados

    resumen_rlg <- objeto |>
      dplyr::group_by(dplyr::across(dplyr::all_of(factores_sin_luz))) |>
      dplyr::summarise(
        RLG_media = mean(.data$RLG, na.rm = TRUE),
        .groups = "drop"
      )

    tabla_fmt <- data.frame(resumen_rlg[factores_sin_luz], check.names = FALSE)
    tabla_fmt$RLG <- vapply(resumen_rlg$RLG_media, function(m) {
      if (is.na(m) || is.nan(m)) return("-")
      return(as.character(round(m, 2)))
    }, character(1))

    tabla_fmt <- .centrar_columna_texto(tabla_fmt, "RLG", "RLG")

    cli::cli_inform(c(
      "v" = "Reporte de fotoblastismo generado.",
      "i" = "RLG: 1 = Fotoblastismo positivo | 0.5 = Indiferente | 0 = Fotoblastismo negativo"
    ))

    return(as.data.frame(tabla_fmt))
  }
}

# ==============================================================================
# FUNCIÓN AUXILIAR INTERNA: Centrado estético de columnas
# ==============================================================================
.centrar_columna_texto <- function(df, col_interna, nom_visible) {
  valores <- as.character(df[[col_interna]])
  valores[is.na(valores) | valores == "NA" | valores == ""] <- "-"

  max_largo <- max(nchar(valores), na.rm = TRUE)
  ancho_celda <- max(max_largo, nchar(nom_visible), 12)

  valores_centrados <- vapply(valores, function(texto) {
    espacios_totales <- ancho_celda - nchar(texto)
    izq <- floor(espacios_totales / 2)
    der <- espacios_totales - izq
    paste0(paste(rep(" ", izq), collapse = ""), texto, paste(rep(" ", der), collapse = ""))
  }, character(1), USE.NAMES = FALSE)

  df[[col_interna]] <- valores_centrados

  espacios_titulo <- ancho_celda - nchar(col_interna)
  izq_tit <- floor(espacios_titulo / 2)
  der_tit <- ancho_celda - (nchar(col_interna) + izq_tit)
  titulo_centrado <- paste0(paste(rep(" ", izq_tit), collapse = ""),
                            col_interna,
                            paste(rep(" ", der_tit), collapse = ""))

  names(df)[names(df) == col_interna] <- titulo_centrado
  return(df)
}
