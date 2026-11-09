# =============================================================================
# 03_metricas_globales.R
#
# Descripción:
#   Calcula métricas globales de estabilidad de la estructura de comunidades
#   entre ventanas consecutivas (t -> t+1) y genera sus gráficas:
#
#     - Métrica propuesta: media armónica de la precisión y el recall del
#       emparejamiento por Jaccard (se descartan las comunidades triviales).
#     - Métrica ajustada: la métrica propuesta multiplicada por la media
#       geométrica de la proporción de vértices en comunidades no triviales en
#       t y en t+1.
#     - Vértices inactivos (singletons estables): proporción de vértices que
#       están en la comunidad 0 en t y siguen ahí en t+1.
#     - Proporción de comunidades no triviales por ventana.
#     - Información mutua normalizada (NMI, implementada en igraph) entre
#       particiones consecutivas.
#
# Entradas:
#   - comunidadesporventana.rds : salida de 01_optimizar_particiones.R
#   - jaccardstotales.rds       : salida de 01_optimizar_particiones.R
#   - recalls.rds               : salida de 01_optimizar_particiones.R
#   - Fechas.csv                : salida de 01_optimizar_particiones.R
#                                 (fechas en formato AAAA-MM-DD)
#   - informacionTrayectorias.rds : salida de 02_seguimiento_completo.R
#   - lista de vectores de membresía : la misma entrada de 01 (para el NMI)
#
# Salidas:
#   - figuras/*.png (si `guardar_figuras` es TRUE)
#
# Uso:
#   1. Ajustar la sección "Configuración".
#   2. Ejecutar:  source("03_metricas_globales.R")
#
# Autor: Diego Antonio Rodríguez Sánchez
# Licencia: MIT (ver LICENSE)
# =============================================================================

# ------------------------------- Configuración -------------------------------
library(ggplot2)
library(tidyr)

archivo_ventanas     <- "comunidadesporventana.rds"        # Salida de 01
archivo_jaccards     <- "jaccardstotales.rds"              # Salida de 01
archivo_recalls      <- "recalls.rds"                      # Salida de 01
archivo_fechas       <- "Fechas.csv"                       # Salida de 01
archivo_trayectorias <- "informacionTrayectorias.rds"      # Salida de 02
archivo_membresias   <- "lista_de_Vectores_de_membresia.rds"  # Misma entrada de 01
indice_lista         <- 1                                  # Posición dentro del RDS (igual que en 01)

umbral_trivial <- 3   # Una comunidad con este número de vértices o menos es trivial

# Fechas que se marcan con una línea vertical en las gráficas (opcional).
# Ejemplo: as.Date(c("2016-01-15", "2019-10-15")). Vacío: no se marca nada.
fechas_resaltar <- as.Date(character(0))
# Fechas de corte entre periodos para la gráfica de métrica e inactivos (opcional).
fechas_periodos <- as.Date(character(0))

guardar_figuras <- TRUE   # TRUE: guarda cada gráfica como PNG
carpeta_figuras <- "figuras"

# ------------------------------- Carga de datos ------------------------------
for (f in c(archivo_ventanas, archivo_jaccards, archivo_recalls, archivo_fechas,
            archivo_trayectorias, archivo_membresias)) {
  if (!file.exists(f)) stop("No se encontró el archivo de entrada: ", f)
}
ventanas     <- readRDS(archivo_ventanas)
jaccards     <- readRDS(archivo_jaccards)
recalls      <- readRDS(archivo_recalls)
trayectorias <- readRDS(archivo_trayectorias)
miLista      <- readRDS(archivo_membresias)
miLista      <- miLista[[indice_lista]]   # Depende del formato del RDS
dffechas     <- read.csv(archivo_fechas, stringsAsFactors = FALSE)

# ================================== Funciones ================================
#' Filtra los cálculos de Jaccard best match y de recall
#'
#' Elimina los cálculos sobre comunidades triviales (de n vértices o menos).
#'
#' @param jaccards Lista `jaccardstotales.rds` (o `recalls.rds`).
#' @param ventanas Lista `comunidadesporventana.rds`.
#' @param tof      TRUE filtra para jaccards, FALSE filtra para recalls.
#' @param n        Cantidad máxima de vértices de una comunidad trivial.
#' @return La lista sin los cálculos descartados.
filtroSingletons <- function(jaccards, ventanas, tof, n) {
  auxiliar <- jaccards
  if (tof) {
    for (i in 1:length(jaccards)) {
      for (j in 1:length(jaccards[[i]])) {
        x <- jaccards[[i]][[j]][2]
        # OJO: `length(... <= n)` cuenta los elementos de la comparación, no los
        # compara con n; la segunda condición siempre es verdadera. Revisar.
        if (length(ventanas[[i]][[j]]) <= n && length(ventanas[[i + 1]][[x]] <= n)) {
          auxiliar[[i]][[j]] <- NA
        }
      }
      auxiliar[[i]] <- auxiliar[[i]][!is.na(auxiliar[[i]])]
    }
    return(auxiliar)
  } else {
    for (i in 1:length(jaccards)) {
      for (j in 1:length(jaccards[[i]])) {
        x <- jaccards[[i]][[j]][2]
        # OJO: mismo caso que arriba.
        if (length(ventanas[[i + 1]][[j]]) <= n && length(ventanas[[i]][[x]] <= n)) {
          auxiliar[[i]][[j]] <- NA
        }
      }
      auxiliar[[i]] <- auxiliar[[i]][!is.na(auxiliar[[i]])]
    }
    return(auxiliar)
  }
}

#' Métrica propuesta (descarta singletons)
#'
#' @param jaccardstotales Lista `jaccardstotales.rds` filtrada.
#' @param recalljaccards  Lista `recalls.rds` filtrada.
#' @return Vector con el valor de la métrica para cada par de ventanas t, t+1.
metricaPropuesta <- function(jaccardstotales, recalljaccards) {
  vector_fac <- c()
  for (i in 1:length(jaccardstotales)) {
    c <- 0
    for (j in 1:length(jaccardstotales[[i]])) {
      c <- c + jaccardstotales[[i]][[j]][3]
    }
    precision <- c / length(jaccardstotales[[i]])
    c2 <- 0
    for (j in 1:length(recalljaccards[[i]])) {
      c2 <- c2 + recalljaccards[[i]][[j]][3]
    }
    recall <- c2 / length(recalljaccards[[i]])
    fac <- (2 * precision * recall) / (precision + recall)
    vector_fac <- c(vector_fac, fac)
  }
  return(vector_fac)
}

#' Métrica ajustada
#'
#' Considera la media geométrica entre la proporción de vértices en
#' comunidades no triviales en t y en t+1.
#'
#' @param ventanas Lista `comunidadesporventana.rds`.
#' @param metrica  Vector de la métrica propuesta (`metricaPropuesta`).
#' @param n        Orden de la red (número de vértices).
#' @param k        Número de vértices para descartar una comunidad.
#' @return Vector con el valor de la métrica ajustada para cada par t, t+1.
metricaAjustadaPropuesta <- function(ventanas, metrica, n, k) {
  proporcionNoSolitarios <- list()
  for (i in 1:(length(ventanas))) {
    c <- 0
    for (j in 1:length(ventanas[[i]])) {
      if (length(ventanas[[i]][[j]]) <= k) {
        c <- c + length(ventanas[[i]][[j]])
      }
    }
    c <- (n - c) / n
    proporcionNoSolitarios[[i]] <- c
  }
  mediageometrica <- c()
  for (i in 1:(length(ventanas) - 1)) {
    c <- sqrt(proporcionNoSolitarios[[i]] * proporcionNoSolitarios[[i + 1]])
    mediageometrica <- c(mediageometrica, c)
  }
  metricamodificada <- c()
  for (i in 1:(length(ventanas) - 1)) {
    c <- metrica[[i]] * mediageometrica[i]
    metricamodificada <- c(metricamodificada, c)
  }
  return(metricamodificada)
}

#' Proporción de vértices que son comunidad 0 y lo siguen siendo en t+1
#'
#' @param comunidad_0 Primera trayectoria de `informacionTrayectorias.rds`.
#' @param n           Número de vértices por ventana.
#' @return Vector con la proporción para cada par t, t+1.
singletonsEstables <- function(comunidad_0, n) {
  solitarios <- c()
  for (i in 1:(length(comunidad_0) - 1)) {
    auxiliar <- comunidad_0[[i]][[3]]
    auxiliar2 <- comunidad_0[[i + 1]][[3]]
    c <- length(intersect(auxiliar, auxiliar2)) / n
    solitarios <- c(solitarios, c)
  }
  return(solitarios)
}

#' Proporción de comunidades no triviales
#'
#' @param ventanas Lista `comunidadesporventana.rds`.
#' @param k        Cantidad máxima de vértices de una comunidad trivial.
#' @return Vector con la proporción para cada ventana.
comunidadNoTrivial <- function(ventanas, k) {
  proporcion <- c()
  for (i in 1:length(ventanas)) {
    c <- 0
    n <- length(ventanas[[i]])
    for (j in 1:n) {
      if (length(ventanas[[i]][[j]]) > k) {
        c <- c + 1
      }
    }
    c <- c / n
    proporcion <- c(proporcion, c)
  }
  return(proporcion)
}

#' Información mutua normalizada entre particiones consecutivas
#'
#' Usa la implementación de igraph sobre las membresías originales. Supone que
#' los vértices están en el mismo orden en todas las ventanas.
#'
#' @param miLista Lista de vectores de membresía, uno por ventana.
#' @return Vector de longitud (número de ventanas - 1).
nmiEntreVentanas <- function(miLista) {
  n <- length(miLista)
  nmi_values <- numeric(n - 1)
  for (i in seq_len(n - 1)) {
    m1 <- miLista[[i]]
    m2 <- miLista[[i + 1]]
    naux <- min(length(m1), length(m2))
    nmi_values[i] <- igraph::compare(
      m1[seq_len(naux)],
      m2[seq_len(naux)],
      method = "nmi"
    )
  }
  return(nmi_values)
}

# ================================== Gráficas =================================
#' Muestra una gráfica y, si `guardar_figuras` es TRUE, la guarda como PNG
#' @param grafica Objeto ggplot.
#' @param nombre  Nombre del archivo, sin extensión.
muestra_grafica <- function(grafica, nombre, ancho = 10, alto = 5) {
  print(grafica)
  if (guardar_figuras) {
    dir.create(carpeta_figuras, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(carpeta_figuras, paste0(nombre, ".png")),
           plot = grafica, width = ancho, height = alto, dpi = 150)
  }
  invisible(grafica)
}

#' Tema común de las gráficas con fechas en el eje x
tema_fechas <- function() {
  theme_bw() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.minor.x = element_line(color = "gray60", linewidth = 0.4),
      panel.grid.major.x = element_line(color = "gray40", linewidth = 0.4)
    )
}

#' Líneas verticales (y etiquetas) en las fechas marcadas
#'
#' @param fechas_marca Fechas (Date) que se marcan. Si está vacío, no agrega nada.
#' @param y            Altura de las etiquetas.
#' @param formato      Formato de la etiqueta (ver `format.Date`).
#' @param angulo,vjust,tamano Apariencia de las etiquetas.
#' @param etiquetas    TRUE agrega las etiquetas con la fecha.
marcas_fechas <- function(fechas_marca, y, formato = "%d %b %Y", angulo = 0,
                          vjust = 0.9, tamano = 3.5, etiquetas = TRUE) {
  if (length(fechas_marca) == 0) return(NULL)
  capas <- list(geom_vline(xintercept = as.numeric(fechas_marca),
                           color = "black", linewidth = 0.8))
  if (etiquetas) {
    capas <- c(capas, list(
      annotate("text", x = fechas_marca, y = y,
               label = format(fechas_marca, formato),
               angle = angulo, hjust = 0, vjust = vjust,
               size = tamano, color = "black")))
  }
  capas
}

#' Gráfica de una serie en el tiempo
#' @param df             data.frame con las columnas `tiempo` y `valor`.
#' @param color          Color de la línea.
#' @param etiqueta_y     Texto del eje y.
#' @param fechas_grafica Fechas (Date) de las marcas principales del eje x.
#' @param titulo,etiqueta_x Título y texto del eje x (opcionales).
#' @param marcas         Capas de `marcas_fechas` (opcional).
grafica_serie <- function(df, color, etiqueta_y, fechas_grafica,
                          titulo = NULL, etiqueta_x = NULL, marcas = NULL) {
  ggplot(df, aes(x = tiempo, y = valor)) +
    geom_line(color = color, linewidth = 0.6) +
    geom_point(color = "black", size = 1) +
    marcas +
    scale_x_date(breaks = fechas_grafica, minor_breaks = df$tiempo,
                 date_labels = "%d %b %Y") +
    labs(title = titulo, x = etiqueta_x, y = etiqueta_y) +
    tema_fechas()
}

#' Gráfica comparada de varias series con valores entre 0 y 1
#' @param df_ancho       data.frame con `tiempo` y una columna por serie.
#' @param colores        Vector con nombre: color de cada serie.
#' @param etiquetas      Vector con nombre: texto de la leyenda de cada serie.
#' @param fechas_grafica Fechas (Date) de las marcas principales del eje x.
#' @param marcas         Capas de `marcas_fechas` (opcional).
#' @param referencias    Capas adicionales, como líneas horizontales (opcional).
grafica_comparada <- function(df_ancho, colores, etiquetas, fechas_grafica,
                              marcas = NULL, referencias = NULL) {
  df_largo <- pivot_longer(df_ancho, cols = -tiempo,
                           names_to = "variable", values_to = "valor")
  ggplot(df_largo, aes(x = tiempo, y = valor, color = variable)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1) +
    referencias +
    marcas +
    scale_x_date(breaks = fechas_grafica, minor_breaks = df_largo$tiempo,
                 date_labels = "%d %b %Y") +
    scale_color_manual(values = colores, labels = etiquetas) +
    labs(x = NULL, y = "Valor", color = "Métrica") +
    ylim(0, 1) +
    tema_fechas()
}

# ================================== Cálculos =================================
numeroVentanas <- length(ventanas)       # Número de ventanas
numeroVertices <- length(miLista[[1]])   # Número de vértices

if (nrow(dffechas) < numeroVentanas) {
  stop(archivo_fechas, " tiene menos fechas que ventanas (", numeroVentanas, ").")
}
fechas_todas <- as.Date(dffechas[1:numeroVentanas, 2])
if (anyNA(fechas_todas)) {
  stop("Las fechas de ", archivo_fechas, " deben tener formato AAAA-MM-DD.")
}
# Las métricas comparan t con t+1: se asignan a la ventana t+1
fechas        <- fechas_todas[2:numeroVentanas]
fechasGrafica <- fechas[seq(1, length(fechas), by = 4)]

# Preparamos los datos
jaccardsFiltro <- filtroSingletons(jaccards, ventanas, TRUE,  umbral_trivial)
recallsFiltro  <- filtroSingletons(recalls,  ventanas, FALSE, umbral_trivial)

# Métrica propuesta y métrica ajustada
metrica         <- metricaPropuesta(jaccardsFiltro, recallsFiltro)
metricaAjustada <- metricaAjustadaPropuesta(ventanas, metrica, numeroVertices, umbral_trivial)

# NMI, vértices inactivos y comunidades no triviales
nmi_values <- nmiEntreVentanas(miLista)
singEsta   <- singletonsEstables(trayectorias[[1]], numeroVertices)
propcom    <- comunidadNoTrivial(ventanas, umbral_trivial)

if (length(singEsta) != length(fechas) || length(nmi_values) != length(fechas)) {
  stop("Las series no tienen la misma longitud: revisa que la comunidad 0 de ",
       archivo_trayectorias, " tenga una entrada por ventana.")
}

# ============================ Gráficas por métrica ===========================
# Métrica ajustada
datafmetAj <- data.frame(tiempo = fechas, valor = metricaAjustada)
muestra_grafica(
  grafica_serie(datafmetAj, "red", "Métrica ajustada", fechasGrafica,
                marcas = marcas_fechas(fechas_resaltar, max(datafmetAj$valor, na.rm = TRUE),
                                       angulo = -90)),
  "metrica_ajustada")

# Métrica propuesta (sin ajuste)
datafmet <- data.frame(tiempo = fechas, valor = metrica)
muestra_grafica(
  grafica_serie(datafmet, "blue", "Métrica propuesta", fechasGrafica,
                marcas = marcas_fechas(fechas_resaltar, 0, etiquetas = FALSE)),
  "metrica_propuesta")

# NMI entre particiones consecutivas
df_nmi <- data.frame(tiempo = fechas, valor = nmi_values)
muestra_grafica(
  grafica_serie(df_nmi, "green", "NMI (1 = idénticas, 0 = independientes)", fechasGrafica,
                titulo = "NMI entre particiones consecutivas",
                etiqueta_x = "Ventana t → t+1",
                marcas = marcas_fechas(fechas_resaltar, max(df_nmi$valor, na.rm = TRUE),
                                       formato = " %b %Y")),
  "nmi")

# Vértices inactivos (singletons estables)
df_singEsta <- data.frame(tiempo = fechas, valor = singEsta)
muestra_grafica(
  grafica_serie(df_singEsta, "orange", "Proporción", fechasGrafica,
                titulo = "Proporción de vértices inactivos",
                etiqueta_x = "Comparación t → t+1",
                marcas = marcas_fechas(fechas_resaltar, max(df_singEsta$valor, na.rm = TRUE),
                                       formato = " %b %Y")),
  "vertices_inactivos")

# Comunidades no triviales (una por cada ventana, no por cada par)
df_propcom <- data.frame(tiempo = fechas_todas, valor = propcom)
muestra_grafica(
  grafica_serie(df_propcom, "purple", "Proporción",
                fechas_todas[seq(1, length(fechas_todas), by = 4)],
                titulo = "Proporción de comunidades no triviales",
                etiqueta_x = "Ventana t"),
  "comunidades_no_triviales")

# ============================ Gráficas comparadas ============================
df_ancho <- data.frame(
  tiempo          = fechas,
  MetricaAjustada = metricaAjustada,
  Inactivos       = singEsta,
  NMI             = nmi_values
)
muestra_grafica(
  grafica_comparada(
    df_ancho,
    colores   = c(MetricaAjustada = "red", Inactivos = "orange", NMI = "green"),
    etiquetas = c(MetricaAjustada = "Métrica propuesta ajustada",
                  Inactivos = "Vértices inactivos",
                  NMI = "Información mutua normalizada"),
    fechasGrafica,
    marcas = marcas_fechas(fechas_resaltar, max(unlist(df_ancho[-1]), na.rm = TRUE),
                           formato = " %b %Y", vjust = 0, tamano = 3.3)),
  "metricas_comparadas")

# Métrica e inactivos con medias y medianas (útil para definir periodos en 05)
df_ancho2 <- df_ancho[, c("tiempo", "MetricaAjustada", "Inactivos")]
mediamet    <- mean(metricaAjustada,   na.rm = TRUE)
medianamet  <- median(metricaAjustada, na.rm = TRUE)
mediasing   <- mean(singEsta,          na.rm = TRUE)
medianasing <- median(singEsta,        na.rm = TRUE)
message(sprintf("Métrica ajustada: media %.3f, mediana %.3f | Inactivos: media %.3f, mediana %.3f",
                mediamet, medianamet, mediasing, medianasing))
referencias <- list(
  geom_hline(yintercept = mediamet,    color = "darkred",    linetype = "dashed", linewidth = 0.8),
  geom_hline(yintercept = medianamet,  color = "pink",       linetype = "dotted", linewidth = 0.8),
  geom_hline(yintercept = mediasing,   color = "darkorange", linetype = "dashed", linewidth = 0.8),
  geom_hline(yintercept = medianasing, color = "coral",      linetype = "dotted", linewidth = 0.8)
)
muestra_grafica(
  grafica_comparada(
    df_ancho2,
    colores   = c(MetricaAjustada = "red", Inactivos = "orange"),
    etiquetas = c(MetricaAjustada = "Métrica propuesta ajustada",
                  Inactivos = "Vértices inactivos"),
    fechasGrafica,
    marcas = marcas_fechas(fechas_periodos, max(unlist(df_ancho2[-1]), na.rm = TRUE),
                           formato = " %b %Y", vjust = 0, tamano = 3.3),
    referencias = referencias),
  "metrica_inactivos_periodos")
