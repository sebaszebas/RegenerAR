# ==============================================================================
# DECLARACIÓN DE VARIABLES GLOBALES
# Silencia las NOTEs de R CMD check sobre variables usadas en NSE (dplyr, etc.)
# ==============================================================================
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # Pronombres y operadores de tidy evaluation
    ".data", ":=",

    # Variables temporales de post_hoc_germinacion
    "Grupo", "media_cruda", "sd_cruda", "n_rep", "media_para_orden",

    # Variables de preparar_datos_germinacion
    "Exitos", "Fracasos", "Total_Germinadas", "PG",
    "t0", "t25", "t50", "t75", "RLG",

    # Variables de ajustar_modelo
    "nombres_raw", "filas_validas", "n_replicas",

    # Variables de reportar_resultados
    "TMG", "IVG", "TMG_let", "IVG_let",
    "t0_m", "t0_s", "t25_m", "t25_s", "t50_m", "t50_s", "t75_m", "t75_s",
    "RLG_media",

    # Funciones base que R a veces confunde
    "head", "sd", "n"
  ))
}
