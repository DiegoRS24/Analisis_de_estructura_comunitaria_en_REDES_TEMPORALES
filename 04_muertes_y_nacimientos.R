# =============================================================================
# 04_muertes_y_nacimientos.R
#
# Descripción:
#   Detecta y clasifica los eventos de nacimiento y muerte de las comunidades
#   temporales de una red por ventanas, y calcula sus tasas de variación.
#   Realiza tres tareas:
#
#     1. Detección de eventos: genera la información de los eventos de
#        nacimiento y muerte de las comunidades temporales (método GED, con
#        umbrales alpha y beta). El resultado se guarda en `eventos.rds` para
#        no recalcularlo cada vez (ver `recalcular_eventos`).
#     2. Estadísticas de eventos: conteos y gráficas de muertes y nacimientos,
#        totales y por tipo de evento.
#     3. Tasa de variación: tasas de variación de muertes y nacimientos
#        (un CSV y una gráfica por serie).
#
# Tipos de evento:
#   Muerte (a):      -1 sin muerte registrada | 0 absorción | 1 división |
#                     2 fragmentación | 3 Teseo
#   Nacimiento (b):  -1 en la ventana 1 | 0 división con muerte del origen |
#                     1 división con supervivencia | 2 fusión |
#                     3 formación difusa | 4 Teseo
#
# Entradas:
#   - informacionTrayectorias.rds : salida de 02_seguimiento_completo.R
#   - comunidadesporventana.rds   : salida de 01_optimizar_particiones.R
#   - Fechas.csv                  : salida de 01_optimizar_particiones.R
#                                   (fechas en formato AAAA-MM-DD)
#
# Salidas:
#   - eventos.rds
#   - TasaInterAnualMyNs.csv
#   - TasaInterQuincenalMyNs.csv
#   - figuras/*.png (si `guardar_figuras` es TRUE)
#
# Uso:
#   1. Ajustar la sección "Configuración".
#   2. Ejecutar:  source("04_muertes_y_nacimientos.R")
#
# Autor: Diego Antonio Rodríguez Sánchez
# Licencia: MIT (ver LICENSE)
# =============================================================================

# ------------------------------- Configuración -------------------------------
library(ggplot2)
library(tidyr)

archivo_trayectorias <- "informacionTrayectorias.rds"  # Salida de 02
archivo_ventanas     <- "comunidadesporventana.rds"    # Salida de 01
archivo_fechas       <- "Fechas.csv"                   # Salida de 01
archivo_eventos      <- "eventos.rds"                  # Salida de este script

recalcular_eventos <- TRUE   # FALSE: lee `archivo_eventos` y omite la parte 1

umbral_jaccard   <- 0.65     # Debe coincidir con el umbral usado en el script 02
alpha_muerte     <- 0.7      # Umbral alpha de GED para muertes (x1)
beta_muerte      <- 0.8      # Umbral beta de GED para muertes (x1 + x2)
alpha_nacimiento <- 0.7      # Umbral alpha de GED para nacimientos (y1)
beta_nacimiento  <- 0.8      # Umbral beta de GED para nacimientos (y1 + y2)

periodo_anual <- 24          # Ventanas que equivalen a un año (24 si son quincenas)
periodo_corto <- 1           # Ventanas a retroceder para la tasa de corto plazo
archivo_tasa_anual <- "TasaInterAnualMyNs.csv"
archivo_tasa_corta <- "TasaInterQuincenalMyNs.csv"

guardar_figuras <- TRUE      # TRUE: guarda cada gráfica como PNG
carpeta_figuras <- "figuras"

# ------------------------------- Carga de datos ------------------------------
for (f in c(archivo_ventanas, archivo_fechas)) {
  if (!file.exists(f)) stop("No se encontró el archivo de entrada: ", f)
}
ventanas       <- readRDS(archivo_ventanas)
numeroVentanas <- length(ventanas)

dffechas <- read.csv(archivo_fechas, stringsAsFactors = FALSE)
if (nrow(dffechas) < numeroVentanas) {
  stop(archivo_fechas, " tiene menos fechas que ventanas (", numeroVentanas, ").")
}
fechas <- as.Date(dffechas[1:numeroVentanas, 2])
if (anyNA(fechas)) {
  stop("Las fechas de ", archivo_fechas, " deben tener formato AAAA-MM-DD.")
}

# =========================== 1. Detección de eventos =========================
#' Índices de una comunidad de la ventana t con las comunidades de otra ventana
#'
#' @param id           Posición de la comunidad dentro de `ventanat`.
#' @param ventanat     Lista de vectores (comunidades de la ventana de origen).
#' @param ventanatmas1 Lista de vectores (comunidades de la ventana de destino).
#' @param indice       TRUE: índice de Jaccard. FALSE: índice de inclusión
#'                     (|A ∩ B| / |A|, con A la comunidad de origen).
#' @return Lista de vectores `c(id, j, valor)` solo con valores mayores que 0.
indicesComunidad <- function(id, ventanat, ventanatmas1, indice) {
  indices <- list()
  y <- 1
  for (j in 1:length(ventanatmas1)) {
    if (indice) {
      x <- length(intersect(ventanat[[id]], ventanatmas1[[j]])) /
           length(union(ventanat[[id]],     ventanatmas1[[j]]))
    } else {
      x <- length(intersect(ventanat[[id]], ventanatmas1[[j]])) /
           length(ventanat[[id]])
    }
    if (x > 0) {
      indices[[y]] <- c(id, j, x)
      y <- y + 1
    }
  }
  indices
}

#' Ventanas de vida de cada comunidad temporal
#'
#' @param trayectorias Lista de trayectorias (`informacionTrayectorias.rds`).
#' @return Lista de vectores `c(x, y, z)`: ventana de la primera aparición,
#'         ventana de la última aparición e ID de la comunidad temporal.
informacionComunidad <- function(trayectorias) {
  informacion <- list()
  contador <- 1
  for (i in 2:length(trayectorias)) {
    aux <- trayectorias[[i]][[1]][[2]]
    c <- c(aux)
    for (j in 1:length(trayectorias[[i]])) {
      if (aux == trayectorias[[i]][[j]][[2]]) {
        aux <- aux + 1
      } else {
        aux2 <- aux - 1
        c <- c(c, aux2, i - 1)
        informacion[[contador]] <- c
        contador <- contador + 1
        aux <- trayectorias[[i]][[j]][[2]]
        c <- c(aux)
        aux <- aux + 1
      }
      if (j == length(trayectorias[[i]])) {
        x <- trayectorias[[i]][[j]][[2]]
        c <- c(c, x, i - 1)
        informacion[[contador]] <- c
        contador <- contador + 1
      }
    }
  }
  informacion
}

#' Encuentra los IDs originales de las comunidades en su primera y última ventana
#'
#' @param ventanas     Comunidades por ventana (salida de 01).
#' @param informacion  Lista obtenida con `informacionComunidad`.
#' @param trayectorias Lista de trayectorias (`informacionTrayectorias.rds`).
#' @return Lista de vectores `c(x, y, z, w)`: ventana de la primera aparición,
#'         ID de la comunidad en esa ventana, ventana de la última aparición e
#'         ID de la comunidad en esa ventana.
traduccion <- function(ventanas, informacion, trayectorias) {
  traducciones <- list()
  for (i in 1:length(informacion)) {
    v <- informacion[[i]]
    comunidad <- v[3] + 1 # Buscamos ID de comunidad al nacer
    ventana <- v[1]
    for (j in 1:length(trayectorias[[comunidad]])) {
      if (ventana == trayectorias[[comunidad]][[j]][[2]]) {
        ventana <- j
        break
      }
    }
    ids <- trayectorias[[comunidad]][[ventana]][[3]]
    nodo <- ids[1] # vértice a buscar
    ventana <- v[1]
    idinicial <- -1
    for (j in 1:length(ventanas[[ventana]])) {
      for (k in 1:length(ventanas[[ventana]][[j]])) {
        if (nodo == ventanas[[ventana]][[j]][k]) {
          idinicial <- j
          break
        }
      }
    }
    comunidad <- v[3] + 1 # Buscamos ID de comunidad al morir
    ventana <- v[2]
    for (j in 1:length(trayectorias[[comunidad]])) {
      if (ventana == trayectorias[[comunidad]][[j]][[2]]) {
        ventana <- j
        break
      }
    }
    ids <- trayectorias[[comunidad]][[ventana]][[3]]
    nodo <- ids[1] # vértice a buscar
    ventana <- v[2]
    idfinal <- -1
    for (j in 1:length(ventanas[[ventana]])) {
      for (k in 1:length(ventanas[[ventana]][[j]])) {
        if (nodo == ventanas[[ventana]][[j]][k]) {
          idfinal <- j
          break
        }
      }
    }
    c <- c(v[1], idinicial, ventana, idfinal)
    traducciones[[i]] <- c
  }
  traducciones
}

#' Calcula x1, x2, y1 e y2 de cada comunidad temporal
#'
#' Auxiliar de la clasificación de eventos. Con índices de inclusión:
#'   x1, x2: los dos mayores índices de la comunidad (en su última aparición)
#'           con las comunidades de la ventana siguiente (muerte).
#'   y1, y2: los dos mayores índices de la comunidad (en su primera aparición)
#'           con las comunidades de la ventana anterior (nacimiento).
#' Son 0 cuando no aplican (sin muerte registrada o nacimiento en la ventana 1).
#'
#' @param ventanas    Comunidades por ventana (salida de 01).
#' @param informacion Lista obtenida con `traduccion`.
#' @return Lista de vectores
#'         `c(ventana_nac, id_nac, ventana_muerte, id_muerte, x1, x2, y1, y2)`.
infoEventos <- function(ventanas, informacion) {
  infoAgregada <- list()
  n <- length(ventanas)
  for (i in 1:length(informacion)) {
    v <- informacion[[i]]
    if (v[3] != n) { # Tipo de muerte (x1 y x2)
      indices <- indicesComunidad(v[4], ventanas[[v[3]]], ventanas[[(v[3] + 1)]], FALSE)
      x1 <- 0
      for (j in seq_along(indices)) {
        if (x1 < indices[[j]][3]) {
          x1 <- indices[[j]][3]
        }
      }
      for (j in seq_along(indices)) {
        if (x1 == indices[[j]][3]) {
          indices <- indices[-j]
          break
        }
      }
      v <- c(v, x1)
      x2 <- 0
      # OJO: aquí la condición es `> 1`; en el cálculo de y2 (nacimientos) es
      # `>= 1`. Revisar si la diferencia es intencional.
      if (length(indices) > 1) {
        for (j in seq_along(indices)) {
          if (x2 < indices[[j]][3]) {
            x2 <- indices[[j]][3]
          }
        }
      }
      v <- c(v, x2)
    } else { # x1 y x2 son 0 si no hay muerte registrada
      v <- c(v, 0, 0)
    }
    # Tipo de nacimiento (y1 e y2)
    if (v[1] != 1) { # Nacimientos fuera de la ventana 1
      indicesNAC <- indicesComunidad(v[2], ventanas[[v[1]]], ventanas[[(v[1] - 1)]], FALSE)
      y1 <- 0
      for (j in seq_along(indicesNAC)) {
        if (y1 < indicesNAC[[j]][3]) {
          y1 <- indicesNAC[[j]][3]
        }
      }
      for (j in seq_along(indicesNAC)) {
        if (y1 == indicesNAC[[j]][3]) {
          indicesNAC <- indicesNAC[-j]
          break
        }
      }
      v <- c(v, y1)
      y2 <- 0
      if (length(indicesNAC) >= 1) {
        for (j in seq_along(indicesNAC)) {
          if (y2 < indicesNAC[[j]][3]) {
            y2 <- indicesNAC[[j]][3]
          }
        }
      }
      v <- c(v, y2)
    } else { # y1 e y2 son 0 si nació en la ventana 1
      v <- c(v, 0, 0)
    }
    infoAgregada[[i]] <- v
  }
  infoAgregada
}

#' Determina si una comunidad murió o nació por Teseo
#'
#' Agrega dos entradas al vector de la comunidad:
#'   a (posición 9):  3 si murió por Teseo, -1 si no tiene muerte registrada,
#'                    0 en otro caso.
#'   b (posición 10): 4 si nació por Teseo, -1 si nació en la ventana 1,
#'                    0 en otro caso.
#'
#' @param trayectoria Vector de una comunidad obtenido con `infoEventos`.
#' @param ventanas    Comunidades por ventana (salida de 01).
#' @param umbral      Umbral de Jaccard usado en el seguimiento (script 02).
eventosTeseo <- function(trayectoria, ventanas, umbral) {
  t <- trayectoria
  n <- length(ventanas)
  if (t[3] < n) { # Muerte
    t <- c(t, 0)
    indices <- indicesComunidad(t[4], ventanas[[ t[3] ]], ventanas[[ t[3] + 1 ]], TRUE)
    for (j in seq_along(indices)) {
      if (indices[[j]][3] > umbral) {
        t[9] <- 3
        break
      }
    }
  } else {
    t <- c(t, -1)
  }
  if (t[1] > 1) { # Nacimiento
    t <- c(t, 0)
    indices <- indicesComunidad(t[2], ventanas[[ t[1] ]], ventanas[[ t[1] - 1 ]], TRUE)
    for (j in seq_along(indices)) {
      if (indices[[j]][3] > umbral) {
        t[10] <- 4
        break
      }
    }
  } else {
    t <- c(t, -1)
  }
  t
}

#' Clasifica los tipos de muerte y nacimiento de las comunidades temporales
#'
#' @param informacionAgregada Lista obtenida con `infoEventos`.
#' @param informacion         Lista obtenida con `informacionComunidad`.
#' @param a_1  Umbral alpha para x1 (muertes).
#' @param b_1  Umbral beta para x1 + x2 (muertes).
#' @param ventanas Comunidades por ventana (salida de 01).
#' @param umbral   Umbral de Jaccard usado en el seguimiento (script 02).
#' @param a_2  Umbral alpha para y1 (nacimientos).
#' @param b_2  Umbral beta para y1 + y2 (nacimientos).
#' @return Lista de vectores `c(x, y, z, x1, x2, y1, y2, a, b)`, donde x, y, z
#'         son los de `informacionComunidad`, x1, x2, y1, y2 los de
#'         `infoEventos`, y a y b los tipos de muerte y de nacimiento
#'         (ver la cabecera del script).
deteccionEventos <- function(informacionAgregada, informacion, a_1, b_1,
                             ventanas, umbral, a_2, b_2) {
  eventos <- list()
  for (i in 1:length(informacionAgregada)) {
    v <- informacionAgregada[[i]]
    eventos[[i]] <- eventosTeseo(v, ventanas, umbral)
  }
  for (i in 1:length(eventos)) {
    if (eventos[[i]][9] == 0) { # muertes
      decision <- 2
      if (eventos[[i]][5] >= a_1) {
        decision <- 0
      } else if ((eventos[[i]][5] + eventos[[i]][6]) >= b_1) {
        decision <- 1
      }
      eventos[[i]][9] <- decision
    }
    if (eventos[[i]][10] == 0) { # nacimientos
      decision <- 3
      if (eventos[[i]][7] >= a_2) {
        t <- eventos[[i]]
        indices <- indicesComunidad(t[2], ventanas[[ t[1] ]], ventanas[[ t[1] - 1 ]], FALSE)
        aux <- 0
        for (j in 1:length(indices)) {
          if (indices[[j]][3] == t[7]) {
            aux <- indices[[j]][2]
            break
          }
        }
        bandera <- FALSE
        for (j in 1:length(informacionAgregada)) {
          if (informacionAgregada[[j]][3] == (t[1] - 1) && informacionAgregada[[j]][4] == aux) {
            bandera <- TRUE
          }
        }
        if (bandera) {
          decision <- 0
        } else {
          decision <- 1
        }
      } else if ((eventos[[i]][7] + eventos[[i]][8]) >= b_2) {
        decision <- 2
      }
      eventos[[i]][10] <- decision
    }
  }
  eventos2 <- list()
  for (i in 1:length(eventos)) {
    x <- informacion[[i]]
    y <- eventos[[i]]
    eventos2[[i]] <- c(x, y[5:length(y)])
  }
  eventos2
}

# ---- Uso ----
if (recalcular_eventos || !file.exists(archivo_eventos)) {
  if (!file.exists(archivo_trayectorias)) {
    stop("No se encontró el archivo de entrada: ", archivo_trayectorias)
  }
  trayectorias <- readRDS(archivo_trayectorias)
  datos               <- informacionComunidad(trayectorias)
  informacionCompleta <- traduccion(ventanas, datos, trayectorias)
  info_eventos        <- infoEventos(ventanas, informacionCompleta)
  eventos <- deteccionEventos(info_eventos, datos,
                              alpha_muerte, beta_muerte, ventanas, umbral_jaccard,
                              alpha_nacimiento, beta_nacimiento)
  saveRDS(eventos, file = archivo_eventos)  # Evita recalcular en cada ejecución
  message("Eventos detectados y guardados en ", archivo_eventos)
} else {
  eventos <- readRDS(archivo_eventos)
  message("Eventos leídos de ", archivo_eventos)
}

# ========================== 2. Estadísticas de eventos =======================
#' Tema común de las gráficas con fechas en el eje x
#' @param color_menor Color de la cuadrícula menor del eje x.
tema_fechas <- function(color_menor = "gray80") {
  theme_bw() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.minor.x = element_line(color = color_menor, linewidth = 0.4),
      panel.grid.major.x = element_line(color = "gray40",    linewidth = 0.4)
    )
}

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

#' Gráfica de varias series de conteos en el tiempo
#' @param df_ancho       data.frame con la columna `tiempo` y una columna por serie.
#' @param colores        Vector con nombre: color de cada serie.
#' @param etiquetas      Vector con nombre: texto de la leyenda de cada serie.
#' @param titulo_leyenda Título de la leyenda.
#' @param fechas         Fechas de las ventanas (Date).
grafica_series <- function(df_ancho, colores, etiquetas, titulo_leyenda, fechas) {
  df_largo <- pivot_longer(df_ancho, cols = -tiempo,
                           names_to = "variable", values_to = "valor")
  ggplot(df_largo, aes(x = tiempo, y = valor, color = variable)) +
    geom_line() + geom_point() +
    scale_x_date(breaks = fechas[seq(4, length(fechas), by = 4)],
                 minor_breaks = df_largo$tiempo, date_labels = "%d %b %Y") +
    scale_color_manual(values = colores, labels = etiquetas) +
    labs(x = NULL, y = "Valor", color = titulo_leyenda) +
    tema_fechas("gray60")
}

#' Gráfica de barras con la suma de eventos por tipo
#' @param sumas  Vector con nombre: total de eventos de cada tipo.
#' @param titulo Título de la gráfica.
grafica_barras <- function(sumas, titulo) {
  df <- data.frame(categoria = factor(names(sumas), levels = names(sumas)),
                   suma = as.numeric(sumas))
  ggplot(df, aes(x = categoria, y = suma, fill = categoria)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(suma, 2)), vjust = -0.5, size = 5) +
    labs(title = titulo, x = "Tipo", y = "Suma total") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.position = "none")
}

#' Conteo de un tipo de evento por ventana
#'
#' @param informacion Lista de vectores obtenida con `deteccionEventos`.
#' @param e  TRUE: nacimientos. FALSE: muertes.
#' @param nE Tipo de evento. Muerte {-1,0,1,2,3}, nacimiento {-1,0,1,2,3,4}.
#'           Con -10 se cuentan todos los nacimientos o todas las muertes.
#' @param n  Número de ventanas de la red.
#' @param k  Periodo de conteo (agrupa los conteos de k ventanas).
#' @return Vector donde la i-ésima entrada es el número de eventos de la
#'         i-ésima ventana.
conteoEvento <- function(informacion, e, nE, n, k) {
  v <- integer(n)
  if (e) { # NACIMIENTO
    for (i in 1:length(informacion)) {
      if (informacion[[i]][9] == nE) {
        x <- informacion[[i]][1]
        aux <- as.numeric(v[x]) + 1
        v[x] <- aux
      } else if (nE == -10) {
        x <- informacion[[i]][1]
        aux <- as.numeric(v[x]) + 1
        v[x] <- aux
      }
    }
  } else { # MUERTE
    for (i in 1:length(informacion)) {
      if (informacion[[i]][8] == nE) {
        x <- informacion[[i]][2]
        aux <- as.numeric(v[x]) + 1
        v[x] <- aux
      } else if (nE == -10 && informacion[[i]][8] != -1) {
        x <- informacion[[i]][2]
        aux <- as.numeric(v[x]) + 1
        v[x] <- aux
      }
    }
  }
  if (k > 1) {
    c <- c()
    contador <- 1
    suma <- 0
    for (i in 1:length(v)) {
      if (i == n) {
        suma <- suma + as.numeric(v[i])
        c <- c(c, suma)
      } else if (contador <= k) {
        suma <- suma + as.numeric(v[i])
        contador <- contador + 1
      } else {
        c <- c(c, suma)
        suma <- as.numeric(v[i])
        contador <- 2
      }
    }
    return(c)
  } else {
    return(v)
  }
}

# ---- Uso ----
# Muertes y nacimientos totales
nacimientos <- conteoEvento(eventos, TRUE,  -10, numeroVentanas, 1)
muertes     <- conteoEvento(eventos, FALSE, -10, numeroVentanas, 1)
df_totales  <- data.frame(tiempo = fechas, Muertes = muertes, Nacimientos = nacimientos)
muestra_grafica(
  grafica_series(df_totales,
                 colores   = c(Muertes = "blue", Nacimientos = "red"),
                 etiquetas = c(Muertes = "Muertes totales", Nacimientos = "Nacimientos totales"),
                 titulo_leyenda = "Evento", fechas = fechas),
  "eventos_totales")

# Muertes por tipo
muertePorabsorcion    <- conteoEvento(eventos, FALSE, 0, numeroVentanas, 1)
muertePordivision     <- conteoEvento(eventos, FALSE, 1, numeroVentanas, 1)
muertePorfragmentacion <- conteoEvento(eventos, FALSE, 2, numeroVentanas, 1)
muertePorTeseo        <- conteoEvento(eventos, FALSE, 3, numeroVentanas, 1)
sinregistromuerte     <- conteoEvento(eventos, FALSE, -1, numeroVentanas, 1)
df_muertes_tipo <- data.frame(
  tiempo        = fechas,
  Absorcion     = muertePorabsorcion,
  Division      = muertePordivision,
  Fragmentacion = muertePorfragmentacion
)
muestra_grafica(
  grafica_series(df_muertes_tipo,
                 colores   = c(Absorcion = "blue", Division = "red", Fragmentacion = "green"),
                 etiquetas = c(Absorcion = "Muerte por absorción",
                               Division = "Muerte por división",
                               Fragmentacion = "Muerte por fragmentación"),
                 titulo_leyenda = "Tipo de muerte", fechas = fechas),
  "muertes_por_tipo")

# Nacimientos por tipo
nacimientosDivisionMuerteOrigen <- conteoEvento(eventos, TRUE, 0, numeroVentanas, 1)
nacimientosPordivision          <- conteoEvento(eventos, TRUE, 1, numeroVentanas, 1)
nacimientosPorFusion            <- conteoEvento(eventos, TRUE, 2, numeroVentanas, 1)
nacimientosFormacionDifusa      <- conteoEvento(eventos, TRUE, 3, numeroVentanas, 1)
nacimientosPorTeseo             <- conteoEvento(eventos, TRUE, 4, numeroVentanas, 1)
sinregistronacimiento           <- conteoEvento(eventos, TRUE, -1, numeroVentanas, 1)
df_nacimientos_tipo <- data.frame(
  tiempo         = fechas,
  DivisionMuerte = nacimientosDivisionMuerteOrigen,
  Division       = nacimientosPordivision,
  Fusion         = nacimientosPorFusion,
  Formacion      = nacimientosFormacionDifusa
)
muestra_grafica(
  grafica_series(df_nacimientos_tipo,
                 colores   = c(DivisionMuerte = "blue", Division = "red",
                               Fusion = "green", Formacion = "purple"),
                 etiquetas = c(DivisionMuerte = "Nacimiento por división del origen",
                               Division = "División del origen sin muerte",
                               Fusion = "Nacimiento por fusión",
                               Formacion = "Formación difusa"),
                 titulo_leyenda = "Tipo de nacimiento", fechas = fechas),
  "nacimientos_por_tipo")

# Distribución de muertes y nacimientos por tipo
sumas_muertes <- c(sum(sinregistromuerte), sum(muertePorabsorcion), sum(muertePordivision),
                   sum(muertePorfragmentacion), sum(muertePorTeseo))
names(sumas_muertes) <- c("Sin registro", "Absorción", "División", "Fragmentación", "Teseo")
muestra_grafica(grafica_barras(sumas_muertes, "Distribución por tipo de muerte"),
                "distribucion_muertes")

sumas_nacimientos <- c(sum(sinregistronacimiento), sum(nacimientosDivisionMuerteOrigen),
                       sum(nacimientosPordivision), sum(nacimientosPorFusion),
                       sum(nacimientosFormacionDifusa), sum(nacimientosPorTeseo))
names(sumas_nacimientos) <- c("Ventana 1", "División c/muerte", "División c/sup",
                              "Fusión", "Formación difusa", "Teseo")
muestra_grafica(grafica_barras(sumas_nacimientos, "Distribución por tipo de nacimiento"),
                "distribucion_nacimientos")

# =========================== 3. Tasa de variación ============================
#' Tasa de variación porcentual de una serie
#'
#' Tasa = (V_final - V_inicial) / V_inicial * 100. Como se necesita un valor
#' inicial de referencia, las primeras `p` ventanas no tienen tasa. Devuelve NA
#' cuando V_inicial = 0.
#'
#' @param datos Vector de datos (número de eventos por ventana).
#' @param p     Ventanas a retroceder para la comparación.
#' @return Lista con la tasa de cada ventana desde la p + 1.
tasasVar <- function(datos, p) {
  serie <- list()
  for (i in (p + 1):length(datos)) {
    if (datos[i - p] == 0) {
      serie <- c(serie, NA)
    } else {
      aux <- (datos[i] - datos[i - p]) / datos[i - p]
      aux <- aux * 100
      serie <- c(serie, aux)
    }
  }
  serie
}

#' Completa las tasas con NA en las ventanas que solo sirven de referencia
#'
#' @param n     Número total de ventanas.
#' @param serie Lista obtenida con `tasasVar`.
#' @return Vector de longitud n con las tasas.
tasasCompletas <- function(n, serie) {
  k <- (n - length(serie))
  vector <- rep(NA, k)
  for (i in 1:length(serie)) {
    aux <- serie[[i]]
    vector <- c(vector, aux)
  }
  vector
}

#' Gráfica de una serie de tasas de variación
#' @param fechas     Fechas de las ventanas con tasa (Date).
#' @param tasas      Vector de tasas.
#' @param color      Color de la línea.
#' @param etiqueta_y Texto del eje y.
grafica_tasa <- function(fechas, tasas, color, etiqueta_y) {
  df <- data.frame(tiempo = as.Date(fechas), tasa = tasas)
  ggplot(df, aes(x = tiempo, y = tasa)) +
    geom_line(color = color, linewidth = 0.7, na.rm = TRUE) +
    geom_point(color = "black", size = 1, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_date(breaks = df$tiempo[seq(1, nrow(df), by = 4)],
                 minor_breaks = df$tiempo, date_labels = "%d %b %Y") +
    labs(x = NULL, y = etiqueta_y) +
    tema_fechas()
}

#' Calcula las tasas de variación de muertes y nacimientos, las grafica y
#' escribe la tabla en un CSV
#'
#' @param muertes,nacimientos Conteos por ventana.
#' @param fechas        Fechas de todas las ventanas (Date).
#' @param p             Ventanas a retroceder para la tasa.
#' @param archivo_salida CSV de salida.
#' @param etiqueta      Texto que describe el periodo (por ejemplo "inter-anual").
procesa_tasas <- function(muertes, nacimientos, fechas, p, archivo_salida, etiqueta) {
  n <- length(fechas)
  tasasMuertes     <- tasasVar(muertes, p)
  tasasNacimientos <- tasasVar(nacimientos, p)
  fechas_tasa      <- fechas[(p + 1):n]
  id <- chartr("-", "_", etiqueta)

  muestra_grafica(
    grafica_tasa(fechas_tasa, unlist(tasasMuertes), "red",
                 paste("Tasa de variación de muertes", etiqueta)),
    paste0("tasa_muertes_", id))
  muestra_grafica(
    grafica_tasa(fechas_tasa, unlist(tasasNacimientos), "blue",
                 paste("Tasa de variación de nacimientos", etiqueta)),
    paste0("tasa_nacimientos_", id))

  tabla <- data.frame(
    Ventana            = fechas,
    Nacimientos        = nacimientos,
    TasaVarNacimientos = tasasCompletas(n, tasasNacimientos),
    Muertes            = muertes,
    TasaVarMuertes     = tasasCompletas(n, tasasMuertes)
  )
  write.csv(tabla, archivo_salida, row.names = FALSE)
  invisible(tabla)
}

# ---- Uso ----
procesa_tasas(muertes, nacimientos, fechas, periodo_anual,
              archivo_tasa_anual, "inter-anual")
procesa_tasas(muertes, nacimientos, fechas, periodo_corto,
              archivo_tasa_corta, "inter-quincenal")
