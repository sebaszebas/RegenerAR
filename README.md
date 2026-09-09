# RegenerAR 

**Herramientas Unificadas para el Análisis y Modelado de Germinación**

[![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Descripción
RegenerAR automatiza el análisis estadístico de ensayos de germinación en R. Calcula índices cinéticos, ajusta modelos GLM/LM con diagnósticos automáticos (sobredispersión, separación perfecta) y genera tablas de comparación múltiple con letras de Tukey listas para publicación.

## Instalación
```r
install.packages("devtools")
devtools::install_github("sebaszebas/RegenerAR")
library(RegenerAR)
```

## Uso Rápido

### 1. Preparar datos y Tabla Maestra
```r
# Generar matriz maestra
matriz <- preparar_datos_germinacion(mis_datos, factores = c("Luz", "Temperatura"))

# Tabla de Óptimos (Germinación + TMG + IVG)
tabla_opt <- reportar_resultados(matriz, tipo = "optimo", firth = TRUE, verbose = FALSE)
print(tabla_opt)
```

### 2. Tablas de Tiempos y Fotoblastismo
```r
# Tabla de percentiles (t0, t25, t50, t75)
tabla_tiempos <- reportar_resultados(matriz, tipo = "percentiles")
print(tabla_tiempos)

# Índice de Fotoblastismo (RLG)
tabla_rlg <- reportar_resultados(matriz, tipo = "fotoblastismo")
print(tabla_rlg)
```

### 3. Modelado Estadístico Detallado
```r
# Modelo de Germinación Final
mod_ger <- ajustar_modelo(matriz, c("Luz", "Temperatura"), "germinacion", firth = TRUE)
summary(mod_ger) # Ver coeficientes y diagnósticos

# Modelo de TMG (ajusta Gamma automáticamente)
mod_tmg <- ajustar_modelo(matriz, c("Luz", "Temperatura"), "TMG")
summary(mod_tmg)
```

## Características Destacadas
- **Motor Inteligente:** Excluye tratamientos sin réplicas, modela el resto y reincorpora los excluidos automáticamente.
- **Fallback Descriptivo:** Si el modelo estadístico falla, calcula Media ± DE sin romper el flujo.
- **Rigor Biológico:** Detecta datos acumulados y aplica blindaje de monotonía en percentiles.
- **Resumen Nativo:** Usa `summary(modelo)` para ver coeficientes y diagnósticos estilo R clásico.

## Documentación
Para ver la guía completa de uso, ejecuta en R:
```r
vignette("guia_uso_RegenerAR", package = "RegenerAR")
```

## Autor
Sebastian R. Zeballos  
*Licencia: GPL-3*
