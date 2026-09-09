#' Post-hoc Unificado para Ensayos de Germinación (Formato Científico)
#'
#' @description
#' Realiza comparaciones múltiples post-hoc sobre el objeto generado por
#' \code{ajustar_modelo()}, aplicando back-transformation cuando corresponde
#' (arcoseno, log) y generando letras de agrupamiento en formato científico.
#'
#' @details
#' ## Lógica de ordenamiento de letras
#' Para variables de tiempo (TMG, t0, t25, t50, t75), las letras se asignan
#' ordenando de menor a mayor (el menor tiempo recibe "a"). Para el resto
#' (germinación, IVG, SYN, etc.), se ordena de mayor a menor.
#'
#' ## Back-transformation automático
#' \itemize{
#'   \item Arcoseno: \code{sin(x)^2} (para SYN, Z_INDEX, RLG)
#'   \item Gamma: medias en escala respuesta (no log)
#'   \item Binomial: probabilidades convertidas a porcentajes
#' }
#'
#' ## Manejo de diseños asimétricos
#' Si un factor fue removido por falta de datos (ej: Luz en oscuridad),
#' las letras se asignan solo para la condición con datos y se reinyectan
#' en la tabla final.
#'
#' @param objeto_modelo Objeto S3 de clase \code{"modelo_germinacion_unificado"}
#'   generado por \code{ajustar_modelo()}.
#' @param alpha Nivel de significancia para comparaciones múltiples. Por defecto 0.05.
#'
#' @return Un data frame con:
#' \describe{
#'   \item{Factores}{Columnas de los factores explicativos}
#'   \item{Valor}{Media ± DE en escala original (o "-" si no hay datos)}
#'   \item{*_let}{Columnas con letras de agrupamiento (una por factor o interacción)}
#' }
#'
#' @importFrom emmeans emmeans
#' @importFrom multcomp cld
#' @importFrom dplyr group_by across all_of summarise filter left_join select arrange desc n
#' @importFrom cli cli_abort cli_inform
#' @export
#'
#' @examples
#' \dontrun{
#' modelo <- ajustar_modelo(matriz_maestra, c("Temperatura", "Luz"), "TMG")
#' posthoc <- post_hoc_germinacion(modelo, alpha = 0.05)
#' print(posthoc)
#' }
post_hoc_germinacion <- function(objeto_modelo, alpha = 0.05) {

  # ============================================================================
  # VALIDACIONES INICIALES
  # ============================================================================
  if (!inherits(objeto_modelo, "modelo_germinacion_unificado")) {
    cli::cli_abort("El objeto debe haber sido generado por {.fn ajustar_modelo}.")
  }

  if (!requireNamespace("emmeans", quietly = TRUE)) {
    cli::cli_abort("El paquete {.pkg emmeans} es necesario. Instálalo con {.code install.packages('emmeans')}.")
  }
  if (!requireNamespace("multcomp", quietly = TRUE)) {
    cli::cli_abort("El paquete {.pkg multcomp} es necesario. Instálalo con {.code install.packages('multcomp')}.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    cli::cli_abort("El paquete {.pkg dplyr} es necesario. Instálalo con {.code install.packages('dplyr')}.")
  }

  # ============================================================================
  # SECCIÓN 1: EXTRACCIÓN DE METADATOS
  # ============================================================================
  firth_en_familia  <- grepl("Firth", objeto_modelo$familia, ignore.case = TRUE)
  modelo            <- objeto_modelo$modelo
  factores          <- objeto_modelo$factores
  factores_activos  <- objeto_modelo$factores_activos
  tipo_modelo       <- objeto_modelo$tipo_modelo
  variable_nombre   <- objeto_modelo$variable
  es_binomial       <- objeto_modelo$es_binomial
  es_gamma          <- objeto_modelo$es_gamma
  datos_originales  <- objeto_modelo$data

  # ============================================================================
  # FUNCIÓN AUXILIAR: Ordenamiento de letras según tipo de variable
  # ============================================================================
  .ordenar_letras <- function(df, col_media = "emmean", variable_nombre) {
    df <- as.data.frame(df)
    columna_letras <- intersect(c(".group", ".letters"), names(df))
    df$Grupo <- if (length(columna_letras) == 0) "a" else trimws(df[[columna_letras]])

    opciones_media <- intersect(c(col_media, "response", "prob", "prediction"), names(df))
    media_activa <- if (length(opciones_media) > 0) opciones_media[1] else names(df)[1]

    # Lógica de inversión fisiológica:
    # TMG, t0-t75 → menor tiempo = "a" (decreasing = FALSE)
    # Germinación, IVG, SYN → mayor valor = "a" (decreasing = TRUE)
    es_variable_tiempo <- variable_nombre %in% c("TMG", "t0", "t25", "t50", "t75")
    orden <- order(df[[media_activa]], decreasing = !es_variable_tiempo, na.last = TRUE)

    df_orden <- df[orden, , drop = FALSE]
    letras_originales <- unique(unlist(strsplit(df_orden$Grupo, split = "")))
    letras_originales <- letras_originales[letras_originales != ""]

    nombres <- stats::setNames(letters[seq_along(letras_originales)], letras_originales)
    df$Grupo <- vapply(df$Grupo, function(grupo) {
      caracteres_nuevos <- nombres[unlist(strsplit(grupo, split = ""))]
      paste(caracteres_nuevos[!is.na(caracteres_nuevos)], collapse = "")
    }, character(1))
    return(df)
  }

  # ============================================================================
  # SECCIÓN 2: BACK-TRANSFORMATION Y RESUMEN DESCRIPTIVO
  # ============================================================================

  # 🔧 CORRECCIÓN CRÍTICA: Usar grepl() para comparación robusta
  if (grepl("Arcoseno", objeto_modelo$transformacion, ignore.case = TRUE)) {
    datos_originales[[variable_nombre]] <- sin(datos_originales[[variable_nombre]])^2
  }

  if (es_binomial) {
    variable_analizada <- "PG_virtual"
    datos_originales$PG_virtual <- (datos_originales$Exitos / (datos_originales$Exitos + datos_originales$Fracasos)) * 100
    titulo_valor_col <- "Germinación (%)"
  } else {
    variable_analizada <- variable_nombre
    titulo_valor_col <- variable_nombre
  }

  lista_niveles <- lapply(datos_originales[factores], function(x) {
    if (is.factor(x)) levels(x) else unique(x)
  })
  combinaciones_totales <- do.call(expand.grid, c(lista_niveles, list(stringsAsFactors = FALSE)))

  tabla_base_resumen <- datos_originales |>
    dplyr::group_by(dplyr::across(dplyr::all_of(factores_activos))) |>
    dplyr::summarise(
      media_cruda = mean(.data[[variable_analizada]], na.rm = TRUE),
      sd_cruda    = stats::sd(.data[[variable_analizada]], na.rm = TRUE),
      n_rep       = dplyr::n(),
      .groups     = "drop"
    )

  tabla_completa <- dplyr::left_join(combinaciones_totales, tabla_base_resumen, by = factores_activos)

  tratamientos_con_datos <- tabla_completa |> dplyr::filter(!is.na(media_cruda) & !is.nan(media_cruda))
  tratamientos_vacios    <- tabla_completa |> dplyr::filter(is.na(media_cruda) | is.nan(media_cruda))

  texto_valor <- character(nrow(tratamientos_con_datos))
  for (i in seq_len(nrow(tratamientos_con_datos))) {
    m <- tratamientos_con_datos$media_cruda[i]
    d <- tratamientos_con_datos$sd_cruda[i]
    n <- tratamientos_con_datos$n_rep[i]

    if (is.na(d) || d == 0 || n < 2) {
      texto_valor[i] <- as.character(round(m, if (es_binomial) 1 else 2))
    } else {
      texto_valor[i] <- paste0(round(m, if (es_binomial) 1 else 2), " ± ", round(d, if (es_binomial) 1 else 2))
    }
  }

  tabla_final <- data.frame(
    tratamientos_con_datos[factores],
    Valor = texto_valor,
    media_para_orden = tratamientos_con_datos$media_cruda,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(tabla_final)[names(tabla_final) == "Valor"] <- titulo_valor_col

  # ============================================================================
  # SECCIÓN 3: EJECUCIÓN DE COMPARACIONES MÚLTIPLES (EMMEANS)
  # ============================================================================

  # 🔧 CORRECCIÓN CRÍTICA: Para Gamma, usar "response" (no "link")
  specs_link <- if (firth_en_familia) "link" else "response"

  if (tipo_modelo == "Interaccion") {
    medias <- emmeans::emmeans(modelo, specs = factores_activos, type = specs_link)
    letras <- multcomp::cld(medias, by = NULL, Letters = letters, adjust = "tukey", alpha = alpha)

    tabla_letras <- as.data.frame(letras) |>
      .ordenar_letras(col_media = if (es_binomial) "prob" else "response", variable_nombre) |>
      dplyr::select(dplyr::all_of(factores_activos), Grupo)

    # 🔧 CORRECCIÓN: Usar collapse para crear un solo nombre
    nombre_let <- paste0(paste(factores_activos, collapse = "_"), "_let")
    names(tabla_letras)[names(tabla_letras) == "Grupo"] <- nombre_let
    tabla_final <- tabla_final |> dplyr::left_join(tabla_letras, by = factores_activos)

  } else if (tipo_modelo == "Aditivo") {
    p_factores <- objeto_modelo$p_factores
    factores_significativos <- factores_activos[!is.na(p_factores) & p_factores < alpha]

    if (length(factores_significativos) > 0) {
      for (f in factores_significativos) {
        medias_f <- emmeans::emmeans(modelo, specs = f, type = specs_link)
        letras_f <- multcomp::cld(medias_f, Letters = letters, adjust = "sidak", alpha = alpha)

        tabla_letras_f <- as.data.frame(letras_f) |>
          .ordenar_letras(col_media = if (es_binomial) "prob" else "response", variable_nombre) |>
          dplyr::select(dplyr::all_of(f), Grupo)

        names(tabla_letras_f)[names(tabla_letras_f) == "Grupo"] <- paste0(f, "_let")
        tabla_final <- tabla_final |> dplyr::left_join(tabla_letras_f, by = f)
      }
    }
  } else {
    # Bloque Un_Factor (incluye asimetría)
    medias <- emmeans::emmeans(modelo, specs = factores_activos, type = specs_link)

    letras <- tryCatch({
      multcomp::cld(medias, Letters = letters, adjust = "tukey", alpha = alpha)
    }, error = function(e) {
      tryCatch({
        multcomp::cld(medias, Letters = letters, adjust = "sidak", alpha = alpha)
      }, error = function(e2) {
        multcomp::cld(medias, Letters = letters, adjust = "none", alpha = alpha)
      })
    })

    tabla_letras <- as.data.frame(letras) |>
      .ordenar_letras(col_media = if (es_binomial) "prob" else "response", variable_nombre) |>
      dplyr::select(dplyr::all_of(factores_activos), Grupo)

    # 🔧 CORRECCIÓN: Crear nombre único
    nombre_let <- paste0(paste(factores_activos, collapse = "_"), "_let")
    names(tabla_letras)[names(tabla_letras) == "Grupo"] <- nombre_let

    # Manejo de factor removido (asimetría)
    factor_removido <- setdiff(factores, factores_activos)

    if (length(factor_removido) > 0) {
      nivel_sobreviviente <- unique(objeto_modelo$data[[factor_removido]])
      tabla_letras[[factor_removido]] <- nivel_sobreviviente[1]
      tabla_final <- tabla_final |> dplyr::left_join(tabla_letras, by = factores)
    } else {
      tabla_final <- tabla_final |> dplyr::left_join(tabla_letras, by = factores_activos)
    }
  }

  columnas_let <- names(tabla_final)[grepl("_let", names(tabla_final))]

  # ============================================================================
  # SECCIÓN 4: REINYECCIÓN DE TRATAMIENTOS VACÍOS
  # ============================================================================
  if (nrow(tratamientos_vacios) > 0) {
    filas_muertas <- tratamientos_vacios[factores]
    filas_muertas[[titulo_valor_col]] <- "-"
    filas_muertas$media_para_orden <- -Inf

    if (length(columnas_let) > 0) {
      for (col in columnas_let) filas_muertas[[col]] <- "-"
    } else {
      nombre_defecto <- paste0(paste(factores_activos, collapse = "_"), "_let")
      filas_muertas[[nombre_defecto]] <- "-"
      columnas_let <- nombre_defecto
    }
    filas_muertas <- filas_muertas[, names(tabla_final), drop = FALSE]
    tabla_final <- rbind(tabla_final, filas_muertas)
  }

  tabla_final <- tabla_final |> dplyr::arrange(dplyr::desc(media_para_orden))

  # ============================================================================
  # SECCIÓN 5: BLINDAJE DE CEROS BIOLÓGICOS
  # ============================================================================
  if (length(columnas_let) > 0) {
    tolerancia_cero <- .Machine$double.eps ^ 0.5
    es_cero_numerico <- !is.na(tabla_final$media_para_orden) & abs(tabla_final$media_para_orden) < tolerancia_cero
    for (col in columnas_let) {
      tabla_final[[col]] <- ifelse(es_cero_numerico, "-", tabla_final[[col]])
    }
  }

  # ============================================================================
  # SECCIÓN 6: LIMPIEZA FINAL
  # ============================================================================
  tabla_final <- tabla_final |> dplyr::select(-media_para_orden)

  cli::cli_inform(c(
    "v" = "Post-hoc completado exitosamente.",
    "i" = "Valores representan Media ± DE en escala original."
  ))

  return(as.data.frame(tabla_final))
}
