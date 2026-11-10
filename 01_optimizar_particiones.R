# =============================================================================
# 01_optimizar_particiones.R
#
# Descripción:
#   Script para el análisis de comunidades en redes temporales (ventanas).
#   A partir de una lista de vectores de membresía (un elemento por ventana),
#   realiza tres tareas:
#
#     1. Extracción de datos: exporta nombres de vértices y fechas de cada
#        ventana a archivos CSV.
#     2. Optimización de membresía: transforma la partición de cada ventana
#        (objeto `membership` de igraph) en vectores de IDs de vértices
#        utilizables por el resto de funciones del repositorio.
#     3. Best Match Jaccard: calcula el mejor emparejamiento por índice de
#        Jaccard entre comunidades de ventanas consecutivas y genera
#        `jaccardstotales.rds` y `recalls.rds` para las métricas del
#        análisis global de la red.
#
# Entradas:
#   - lista_de_Vectores_de_membresia.rds : lista de objetos `membership` de igraph.
#     Nota: names(lista_de_Vectores_de_membresia) debe corresponder a
#     las fechas de cada ventana (cada vector de membresía); 
#     names(lista_de_Vectores_de_membresia[[1]]) debe  corresponder a los nombres de 
#     los nodos de la red.
#
# Salidas:
#   - Fechas.csv
#   - NodosNames.csv
#   - comunidadesporventana.rds
#   - jaccardstotales.rds
#   - recalls.rds
#
# Uso:
#   1. Ajustar la sección "Configuración".
#   2. Ejecutar:  source("01_optimizar_particiones.R")
#
# Autor: Diego Antonio Rodríguez Sánchez
# Licencia: MIT (ver LICENSE)
# =============================================================================

# ------------------------------- Configuración -------------------------------
setwd("C:/Users/diego/Downloads/carpetaPrueba")
archivo_entrada <- "membresias_desestacionalizada.rds"  # RDS de entrada
indice_lista    <- 1                                     # Posición dentro del RDS

archivo_fechas  <- "Fechas.csv"
archivo_nodos   <- "NodosNames.csv"

# ------------------------------- Carga de datos ------------------------------
miLista <- readRDS(archivo_entrada)
miLista <- miLista[[indice_lista]]   # Depende del formato del RDS

# ============================ 1. Extracción de datos =========================
# Fechas (una por ventana)
categorias <- names(miLista)
df_fechas  <- data.frame(x = categorias)
lineas <- c(
  '"","x"',
  paste0('"', seq_along(df_fechas$x), '",', df_fechas$x)
)
writeLines(lineas, archivo_fechas)

# Nombres de los vértices
categorias <- names(miLista[[1]])
df_nodos   <- data.frame(x = categorias)
lineas <- c(
  '"","x"',
  paste0('"', seq_along(df_nodos$x), '","', df_nodos$x, '"')
)
writeLines(lineas, archivo_nodos)

# ========================= 2. Optimización de membresía ======================
#' Separa en vectores las comunidades de una ventana
#'
#' @param obj Objeto `membership` de igraph (una ventana).
#' @param df  data.frame con los nombres e IDs de los vértices.
#' @return Lista de vectores; cada vector contiene los IDs de una comunidad.
obtenComunidadesPorVentana <- function(obj, df) {
  comunidadesVentana <- list()
  x   <- length(unique(obj))
  cont <- 1
  for (j in 1:x) {
    comunidad <- c()
    for (i in 1:length(obj)) {
      if (obj[i] == j && !is.na(names(obj[i]))) {
        for (k in 1:length(df[[1]])) {
          if (names(obj[i]) == df[k, 2]) {
            comunidad <- c(comunidad, df[k, 1])
            break
          }
        }
      }
    }
    if (length(comunidad) != 0) {
      comunidadesVentana[[cont]] <- comunidad
      cont <- cont + 1
    }
  }
  comunidadesVentana
}

# ---- Uso ----
df <- read.csv(archivo_nodos, stringsAsFactors = FALSE)
comunidadesporventana <- list()
for (i in seq_along(miLista)) {
  comunidadesporventana[[i]] <- obtenComunidadesPorVentana(miLista[[i]], df)
  message(sprintf("Ventana %d completada", i))
}
saveRDS(comunidadesporventana, file = "comunidadesporventana.rds")

# ========================= 3. Best Match Jaccard =============================
#' Índices de Jaccard entre cada pareja de comunidades de dos ventanas
#'
#' @param ventanat     Lista de vectores (comunidades en t).
#' @param ventanatmas1 Lista de vectores (comunidades en t+1).
#' @return Lista de vectores `c(i, j, jaccard)` con índices no nulos.
jaccardentreVentanas2 <- function(ventanat, ventanatmas1) {
  indices <- list()
  y <- 1
  for (i in 1:length(ventanat)) {
    for (j in 1:length(ventanatmas1)) {
      x <- length(intersect(ventanat[[i]], ventanatmas1[[j]])) /
           length(union(ventanat[[i]],     ventanatmas1[[j]]))
      indices[[y]] <- c(i, j, x)
      y <- y + 1
    }
  }
  indices1 <- list()
  z <- 1
  for (i in 1:length(indices)) {
    if (indices[[i]][3] != 0) {
      indices1[[z]] <- indices[[i]]
      z <- z + 1
    }
  }
  indices1
}

#' Agrupa las comparaciones por comunidad de origen
#' @param indices Lista de índices de Jaccard entre dos ventanas.
#' @return Lista de listas (una por comunidad en t).
separaComunidad <- function(indices) {
  keys <- sapply(indices, function(v) v[1])
  split(indices, keys)
}

#' Selecciona el mejor match (Jaccard máximo) por comunidad
#' @param indicesSeparados Lista generada por `separaComunidad`.
#' @return Lista con un vector `c(i, j, jaccard_max)` por comunidad.
bestMatchFiltro <- function(indicesSeparados) {
  indicesAltos <- list()
  y <- 1
  for (i in 1:length(indicesSeparados)) {
    x <- 0
    for (j in 1:length(indicesSeparados[[i]])) {
      if (x < indicesSeparados[[i]][[j]][3]) {
        x <- indicesSeparados[[i]][[j]][3]
        indicesAltos[[y]] <- indicesSeparados[[i]][[j]]
      }
    }
    y <- y + 1
  }
  indicesAltos
}

#' Best Match Jaccard masivo (t -> t+1) para todas las ventanas
#' @param ventanas Lista de comunidades por ventana.
jaccardbestMatchMasivo <- function(ventanas) {
  jaccardscompletos <- list()
  for (i in 1:(length(ventanas) - 1)) {
    aux  <- jaccardentreVentanas2(ventanas[[i]], ventanas[[i + 1]])
    aux2 <- separaComunidad(aux)
    jaccardscompletos[[i]] <- bestMatchFiltro(aux2)
    message(sprintf("Jaccard v_%d -> v_%d completado", i, i + 1))
  }
  jaccardscompletos
}

#' Recall masivo (t+1 -> t) para todas las ventanas
#' @param ventanas Lista de comunidades por ventana.
jaccardbestMatchMasivoRecall <- function(ventanas) {
  jaccardscompletos <- list()
  for (i in 1:(length(ventanas) - 1)) {
    aux  <- jaccardentreVentanas2(ventanas[[i + 1]], ventanas[[i]])
    aux2 <- separaComunidad(aux)
    jaccardscompletos[[i]] <- bestMatchFiltro(aux2)
    message(sprintf("Recall v_%d -> v_%d completado", i + 1, i))
  }
  jaccardscompletos
}

# ---- Uso ----
jaccardstotales <- jaccardbestMatchMasivo(comunidadesporventana)
saveRDS(jaccardstotales, file = "jaccardstotales.rds")

recalls <- jaccardbestMatchMasivoRecall(comunidadesporventana)
saveRDS(recalls, file = "recalls.rds")
