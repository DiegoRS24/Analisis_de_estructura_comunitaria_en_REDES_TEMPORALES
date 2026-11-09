# =============================================================================
# 02_seguimiento_completo.R
#
# Descripción:
#   A partir de los archivos generados por 01_optimizar_particiones.R, sigue
#   las comunidades a lo largo de las ventanas y construye las comunidades
#   temporales. Realiza tres tareas:
#
#     1. Cálculo de trayectorias: seguimiento de todas las comunidades
#        temporales que surgen en todas las ventanas, con un umbral de Jaccard
#        y un umbral de contención (criterio contra el problema de Teseo).
#     2. Re-etiquetado y mapa aluvial: genera la información para el generador
#        de mapas aluviales de Mapequation (https://www.mapequation.org/alluvial/)
#        como un archivo .tree por ventana, con un re-etiquetado de las
#        comunidades para una mejor visualización.
#     3. Información de comunidades temporales: para cada comunidad temporal
#        guarda, en cada ventana donde vive, los IDs y los nombres de sus
#        vértices. Genera `informacionTrayectorias.rds`.
#
# Entradas:
#   - jaccardstotales.rds       : salida de 01_optimizar_particiones.R
#   - comunidadesporventana.rds : salida de 01_optimizar_particiones.R
#   - Fechas.csv                : salida de 01_optimizar_particiones.R
#   - NodosNames.csv            : salida de 01_optimizar_particiones.R
#
# Salidas:
#   - informacionTrayectorias.rds
#   - Red/<fecha>.tree (un archivo por ventana, para el mapa aluvial)
#
# Nota sobre la indexación:
#   La lista `informacionTrayectorias` empieza por la comunidad 0, que reúne a
#   los vértices que no pertenecen a ninguna comunidad no trivial. Para
#   consultar una comunidad temporal por su ID del mapa aluvial, usa la entrada
#   i + 1:  informacionTrayectorias[[i + 1]]
#
# Uso:
#   1. Ajustar la sección "Configuración".
#   2. Ejecutar:  source("02_seguimiento_completo.R")
#
# Autor: Diego Antonio Rodríguez Sánchez
# Licencia: MIT (ver LICENSE)
# =============================================================================

# ------------------------------- Configuración -------------------------------
archivo_jaccards  <- "jaccardstotales.rds"        # Salida de 01
archivo_ventanas  <- "comunidadesporventana.rds"  # Salida de 01
archivo_fechas    <- "Fechas.csv"                 # Salida de 01
archivo_nodos     <- "NodosNames.csv"             # Salida de 01
archivo_salida    <- "informacionTrayectorias.rds"
carpeta_tree      <- "Red"                        # Carpeta de los archivos .tree

umbral_jaccard    <- 0.65  # Jaccard mínimo para seguir una comunidad entre ventanas
umbral_contencion <- 0.33  # Contención mínima de la comunidad seguida respecto a
                           # su primera aparición (criterio antiteseo)
umbral_trivial    <- 3     # Una comunidad con este número de vértices o menos es trivial

# ------------------------------- Carga de datos ------------------------------
for (f in c(archivo_jaccards, archivo_ventanas, archivo_fechas, archivo_nodos)) {
  if (!file.exists(f)) stop("No se encontró el archivo de entrada: ", f)
}
jaccards  <- readRDS(archivo_jaccards)
ventanas  <- readRDS(archivo_ventanas)
dffechas  <- read.csv(archivo_fechas, stringsAsFactors = FALSE)
dfnombres <- read.csv(archivo_nodos,  stringsAsFactors = FALSE)

# ========================= 1. Cálculo de trayectorias ========================
#' Busca todas las trayectorias de seguimiento de comunidades dado un umbral
#'
#' Una trayectoria es `c(x, y1, y2, ..., yn)` donde x es la primera ventana de
#' aparición de la comunidad y yi es el ID de la comunidad en la ventana
#' x + (i - 1).
#'
#' @param umbral   Valor mínimo de Jaccard aceptado entre dos comunidades.
#' @param listaBMJ Lista `jaccardstotales.rds`, obtenida en el script 01.
#' @param ventanas Lista `comunidadesporventana.rds`, obtenida en el script 01.
#' @param umbral2  Valor mínimo de contención aceptado de la comunidad seguida
#'                 en la ventana t con la comunidad de la ventana inicial
#'                 (antiteseo).
#' @return Lista de trayectorias.
calculatrayectorias <- function(umbral, listaBMJ, ventanas, umbral2) {
  lista <- listaBMJ
  trayectorias <- list()
  for (i in 1:length(lista)) {
    for (j in 1:length(lista[[i]])) {
      if (lista[[i]][[j]][1] != 0 && lista[[i]][[j]][3] > umbral) {
        x <- lista[[i]][[j]][2] # Segundo ID encontrado
        c <- c(i, j, x)
        if (i < length(lista)) {
          for (k in (i + 1):length(lista)) {
            bandera <- TRUE
            for (l in 1:length(lista[[k]])) {
              if (x == lista[[k]][[l]][1] && lista[[k]][[l]][3] > umbral) {
                comunidadOriginal <- ventanas[[i]][[j]]
                comunidadActual <- ventanas[[k]][[l]]
                jaccard <- length(intersect(comunidadOriginal, comunidadActual)) / length(comunidadOriginal)
                if (jaccard >= umbral2) {
                  x <- lista[[k]][[l]][2] # ID actualizado
                  c <- c(c, x)
                  lista[[k]][[l]] <- c(0, 0, 0)
                  bandera <- FALSE
                  break
                } else if (x == lista[[k]][[l]][1]) {
                  lista[[k]][[l]] <- c(0, 0, 0)
                }
              } else if (x == lista[[k]][[l]][1]) {
                lista[[k]][[l]] <- c(0, 0, 0)
              }
            }
            if (bandera) {
              break
            }
          }
        }
        trayectorias <- c(trayectorias, list(c))
      } else { # Aquí se agregan las comunidades que solo sobreviven una ventana
        if (lista[[i]][[j]][1] != 0) {
          c <- c(i, j)
          trayectorias <- c(trayectorias, list(c))
        }
      }
    }
  }
  # Agregamos las comunidades de la última ventana
  aux <- length(ventanas)
  visitadas <- c()
  for (i in 1:length(trayectorias)) { # Buscamos las ya visitadas
    # Ventana de la última aparición de la comunidad
    vfin <- trayectorias[[i]][1] + length(trayectorias[[i]]) - 2
    if (vfin == aux) {
      idfin <- trayectorias[[i]][[length(trayectorias[[i]])]]
      visitadas <- c(visitadas, idfin)
    }
  }
  traultven <- list()
  for (i in 1:length(ventanas[[aux]])) {
    contador <- 0
    for (j in seq_along(visitadas)) {
      if (i != visitadas[j]) {
        contador <- contador + 1
      } else {
        break
      }
    }
    if (contador == length(visitadas)) {
      c <- c(aux, i)
      traultven <- c(traultven, list(c))
    }
  }
  trayectorias2 <- c(trayectorias, traultven)
  return(trayectorias2)
}

#' Marca las trayectorias de singletons y de comunidades triviales
#'
#' Determina si un seguimiento es siempre sobre comunidades de n vértices o
#' menos (singletons si n = 1). Agrega un elemento al final de la trayectoria:
#'   -1 si es un singleton seguido,
#'   -2 si es una comunidad trivial (n vértices o menos).
#'
#' @param trayectorias Lista de trayectorias.
#' @param ventanas     Lista `comunidadesporventana.rds`.
#' @param n            Número de vértices para descartar una comunidad.
#' @return Lista de trayectorias con las marcas.
descarte <- function(trayectorias, ventanas, n) {
  modificaciones <- trayectorias
  for (i in 1:length(trayectorias)) {
    c <- 0
    x <- trayectorias[[i]][1]
    y <- trayectorias[[i]][2]
    bandera <- FALSE
    if (length(ventanas[[x]][[y]]) <= n) { # Siempre singleton
      for (j in 2:length(trayectorias[[i]])) {
        b <- trayectorias[[i]][j]
        if (length(ventanas[[x]][[b]]) == 1) {
          c <- c + 1
        }
        x <- x + 1
      }
      if ((length(trayectorias[[i]]) - 1) == c) {
        v <- modificaciones[[i]]
        v <- c(v, -1)
        modificaciones[[i]] <- v
      } else {
        bandera <- TRUE
      }
    }
    if (bandera) { # Trayectoria trivial
      c <- 0
      x <- trayectorias[[i]][1]
      for (j in 2:length(trayectorias[[i]])) {
        b <- trayectorias[[i]][j]
        if (length(ventanas[[x]][[b]]) <= n) {
          c <- c + 1
        }
        x <- x + 1
      }
      if ((length(trayectorias[[i]]) - 1) == c) {
        v <- modificaciones[[i]]
        v <- c(v, -2)
        modificaciones[[i]] <- v
      }
    }
  }
  return(modificaciones)
}

# ---- Uso ----
trayectorias <- calculatrayectorias(umbral_jaccard, jaccards, ventanas, umbral_contencion)
# Trayectorias con singletons y comunidades triviales marcadas
traysingle <- descarte(trayectorias, ventanas, umbral_trivial)

# ====================== 2. Re-etiquetado y mapa aluvial ======================
#' Relaciona trayectorias que "resurgen"
#'
#' Determina si una trayectoria X comparte similitud con una trayectoria Y con
#' muerte registrada antes del nacimiento de X. Se compara la última comunidad
#' seguida de Y con la primera de X (índice de Jaccard) y la primera de Y con
#' la primera de X (índice de contención). Si se superan ambos umbrales, X se
#' marca con el número de trayectoria de su "ancestro" para asignar el mismo ID
#' de comunidad a ambas.
#'
#' @param tray     Lista de trayectorias con singletons y triviales marcados.
#' @param ventanas Lista `comunidadesporventana.rds`.
#' @param umbral1  Umbral de Jaccard.
#' @param umbral2  Umbral de contención.
#' @return Lista de trayectorias no triviales marcadas con el número de
#'         trayectoria de su ancestro.
resurgimiento <- function(tray, ventanas, umbral1, umbral2) {
  trayectorias <- tray
  for (i in 1:length(trayectorias)) {
    aux <- length(trayectorias[[i]])
    if (trayectorias[[i]][aux] > 0 && i < length(trayectorias)) { # Ni trivial ni singleton
      # Guardamos datos de la comunidad donde nace la trayectoria Y
      py1 <- trayectorias[[i]][1] # ventana
      py2 <- trayectorias[[i]][2] # ID de comunidad
      # Guardamos datos de la comunidad donde muere la trayectoria Y
      uy1 <- py1 + aux - 2 # ventana
      uy2 <- trayectorias[[i]][aux] # ID de comunidad
      for (j in (i + 1):length(trayectorias)) {
        auxJ <- length(trayectorias[[j]])
        if (trayectorias[[i]][1] != trayectorias[[j]][1] && trayectorias[[j]][auxJ] > 0) {
          # Guardamos datos de la comunidad donde nace la trayectoria X
          c <- trayectorias[[j]]
          px1 <- c[1] # ventana
          px2 <- c[2] # ID de la comunidad
          # Calculamos similitud
          jaccard <- length(intersect(ventanas[[px1]][[px2]], ventanas[[uy1]][[uy2]])) / length(union(ventanas[[px1]][[px2]], ventanas[[uy1]][[uy2]]))
          cont <- length(intersect(ventanas[[py1]][[py2]], ventanas[[px1]][[px2]])) / length(ventanas[[py1]][[py2]])
          if (umbral1 < jaccard && umbral2 < cont) {
            rel <- -1 * i
            c <- c(c, rel)
            trayectorias[[j]] <- c
            break
          }
        }
      }
    }
  }
  return(trayectorias)
}

#' Determina a qué comunidad pertenece cada vértice en cada ventana
#'
#' Realiza un re-etiquetado: no conserva los IDs originales del objeto
#' `membership`. Imprime en consola el ID más alto utilizado.
#'
#' @param trayectorias Lista obtenida con `resurgimiento`.
#' @param ventanas     Lista `comunidadesporventana.rds`.
#' @return Lista de vectores `c(x, y, z)`: ID del vértice, número de ventana y
#'         nuevo ID de la comunidad (0 para vértices en comunidades triviales).
infoAluvial <- function(trayectorias, ventanas) {
  # Lista a devolver
  vectores <- list()
  # Contador de vectores agregados
  cont <- 1
  # ID comunidad, suma uno cada que se usa un id; iniciado en 1
  id <- 1
  # Recorremos todas las trayectorias registradas
  listaAuxiliar <- list() # Guardamos cada trayectoria i con su ID asignado
  la <- 1
  for (i in 1:length(trayectorias)) {
    # Determinamos si es una trayectoria trivial y la agregamos a comunidad 0
    if (trayectorias[[i]][length(trayectorias[[i]])] == -1) { # Singletons
      y <- trayectorias[[i]][1] # ventana donde inicia la trayectoria
      z <- trayectorias[[i]][2] # auxiliar para conocer el ID del vértice
      x <- ventanas[[y]][[z]] # ID del vértice
      for (j in 2:(length(trayectorias[[i]]) - 1)) {
        vectores[[cont]] <- c(x, y, 0) # Comunidad de triviales
        y <- y + 1
        cont <- cont + 1
      }
    } else if (trayectorias[[i]][length(trayectorias[[i]])] == -2) { # Trivial
      y <- trayectorias[[i]][1] # ventana donde inicia la trayectoria
      for (j in 2:(length(trayectorias[[i]]) - 1)) {
        x <- trayectorias[[i]][j] # id original de la comunidad
        for (k in 1:length(ventanas[[y]][[x]])) {
          nodo <- ventanas[[y]][[x]][k]
          vectores[[cont]] <- c(nodo, y, 0)
          cont <- cont + 1
        }
        y <- y + 1
      }
    } else { # Trayectorias normales
      if (trayectorias[[i]][length(trayectorias[[i]])] < 0) {
        r <- trayectorias[[i]][length(trayectorias[[i]])] * -1
        auxiliarID <- 0
        for (k in 1:length(listaAuxiliar)) {
          if (listaAuxiliar[[k]][1] == r) {
            auxiliarID <- listaAuxiliar[[k]][2]
            break
          }
        }
        # Trayectorias con ancestro
        y <- trayectorias[[i]][1] # ventana donde inicia la trayectoria
        for (j in 2:(length(trayectorias[[i]]) - 1)) {
          x <- trayectorias[[i]][j] # id original de la comunidad
          for (k in 1:length(ventanas[[y]][[x]])) {
            nodo <- ventanas[[y]][[x]][k]
            vectores[[cont]] <- c(nodo, y, auxiliarID)
            cont <- cont + 1
          }
          y <- y + 1
        }
        c <- c(i, auxiliarID)
        listaAuxiliar[[la]] <- c
        la <- la + 1
      } else {
        # Trayectorias sin ancestro
        y <- trayectorias[[i]][1] # ventana donde inicia la trayectoria
        for (j in 2:length(trayectorias[[i]])) {
          x <- trayectorias[[i]][j] # id original de la comunidad
          for (k in 1:length(ventanas[[y]][[x]])) {
            nodo <- ventanas[[y]][[x]][k]
            vectores[[cont]] <- c(nodo, y, id)
            cont <- cont + 1
          }
          y <- y + 1
        }
        c <- c(i, id)
        listaAuxiliar[[la]] <- c
        la <- la + 1
        id <- id + 1
      }
    }
  }
  print(id - 1)
  return(vectores)
}

#' Escribe un archivo .tree por ventana (formato de Infomap)
#'
#' @param listaVentanas Lista construida en `cluInfomapa`.
#' @param nombresRedes  Fechas correspondientes a cada ventana.
#' @param carpeta       Nombre de la carpeta de salida.
escribirTrees <- function(listaVentanas, nombresRedes, carpeta = "Red") {
  if (!dir.exists(carpeta)) { # Creamos carpeta
    dir.create(carpeta, recursive = TRUE)
    cat("Carpeta creada:", carpeta, "\n")
  }
  for (i in seq_along(listaVentanas)) { # Iteramos cada ventana
    red <- listaVentanas[[i]]
    x <- red[[1]]
    y <- red[[2]]
    z <- red[[3]]
    a <- red[[4]]
    b <- red[[5]]
    # Verificamos tamaños de vectores
    n <- length(x)
    if (length(y) != n || length(z) != n || length(a) != n || length(b) != n) {
      warning(paste("Longitudes diferentes en ventana", i, "-Error"))
      next
    }
    encabezado <- paste0("#Formato .tree de infomap, Vertices ", n) # Escribimos
    lineas <- paste0(
      x, ":", y, " ", round(z, 5), " ", '"', a, '"', " ",
      b
    )
    # Nombre de cada archivo como la fecha de la ventana
    nombre_archivo <- paste0(nombresRedes[[i]], ".tree")
    ruta_completa <- file.path(carpeta, nombre_archivo)
    writeLines(c(encabezado, lineas), con = ruta_completa)
    cat("Guardado:", ruta_completa, "\n")
  }
  cat("\nArchivos .tree creados en:", carpeta, "\n")
}

#' Genera los archivos .tree para el generador de mapas aluviales
#'
#' @param vectores Lista de vectores obtenida con `infoAluvial`.
#' @param nombres  Data frame de `NodosNames.csv`.
#' @param fechas   Data frame de `Fechas.csv`.
#' @param carpeta  Nombre de la carpeta de todos los archivos .tree.
#' @return Lista de listas, una por ventana, con cinco vectores: ID de la
#'         comunidad del vértice, posición del vértice dentro de la comunidad,
#'         peso del vértice (constante), nombre del vértice e ID original del
#'         vértice (definido en `NodosNames.csv`).
cluInfomapa <- function(vectores, nombres, fechas, carpeta = "Red") {
  # Separamos por ventana
  vecOrden <- vectores[order(sapply(vectores, `[`, 2))]
  v <- length(fechas[[1]]) # Número de ventanas
  n <- length(nombres[[1]]) # Número de vértices
  ventanas <- list()
  for (i in 1:v) {
    ventana <- list()
    x <- ((i - 1) * n) + 1
    y <- i * n
    aux <- 1
    for (j in x:y) {
      ventana[[aux]] <- vecOrden[[j]]
      aux <- aux + 1
    }
    final <- ventana[order(sapply(ventana, `[`, 3))]
    ventanas[[i]] <- final
  }
  listaFinal <- list() # Lista que contiene a cada ventana
  for (i in 1:length(ventanas)) {
    id <- ventanas[[i]][[1]][3] # ID de la comunidad
    c <- c() # ID de cada vértice dentro de la comunidad
    c2 <- c() # Comunidad de cada vértice
    c3 <- c() # ID del vértice
    flujo <- rep((1 / n), n)
    contador <- 1
    for (j in 1:length(ventanas[[i]])) {
      if (id == ventanas[[i]][[j]][3]) {
        c <- c(c, contador)
        contador <- contador + 1
        c2 <- c(c2, ventanas[[i]][[j]][3])
        c3 <- c(c3, ventanas[[i]][[j]][1])
      } else {
        contador <- 1
        c <- c(c, contador)
        contador <- contador + 1
        c2 <- c(c2, ventanas[[i]][[j]][3])
        c3 <- c(c3, ventanas[[i]][[j]][1])
        id <- ventanas[[i]][[j]][3]
      }
    }
    names <- nombres[[2]] # Columna que contiene los strings de los nombres
    c4 <- names[c3]
    ven <- list(c2, c, flujo, c4, c3) # Lista de vectores con toda la información
    listaFinal[[i]] <- ven
  }
  escribirTrees(listaFinal, fechas[[2]], carpeta)
  return(listaFinal)
}

# ---- Uso ----
# Relacionamos las trayectorias que coinciden (resurgimientos)
resurg <- resurgimiento(traysingle, ventanas, umbral_jaccard, umbral_contencion)
# Preparamos los datos de todos los vértices
infoNODOS <- infoAluvial(resurg, ventanas)
# Generamos los archivos .tree para subirlos al generador de mapas aluviales.
# La función devuelve también una lista, útil para ggplot pero descartada por la
# cantidad de datos.
aluvial <- cluInfomapa(infoNODOS, dfnombres, dffechas, carpeta_tree)

# ================= 3. Información de comunidades temporales ==================
#' Información de cada comunidad temporal
#'
#' Crea la lista de información de cada comunidad temporal según su ID en el
#' re-etiquetado del mapa aluvial. Cada comunidad es una lista de ventanas, y
#' cada ventana es una lista de 4 datos: fecha, número de ventana, vector de
#' IDs de los vértices de la comunidad en esa ventana y vector de sus nombres.
#'
#' @param vectores          Lista obtenida con `infoAluvial`.
#' @param numeroComunidades ID más alto de las comunidades seguidas.
#' @param nombres           Vector con los nombres de los vértices.
#' @param fechas            Vector con las fechas de las ventanas.
comunidadTrayectoria <- function(vectores, numeroComunidades, nombres, fechas) {
  # separamos por trayectoria/comunidad
  vecOrden <- vectores[order(sapply(vectores, `[`, 3))] # Orden por comunidad
  x <- numeroComunidades
  trayectorias <- list()
  aux <- 1
  for (i in 0:x) {
    listaVectores <- list()
    lugar <- 1
    for (j in aux:length(vectores)) {
      if (vecOrden[[j]][3] == i) {
        listaVectores[[lugar]] <- vecOrden[[j]]
        lugar <- lugar + 1
      } else {
        aux <- j
        break
      }
    }
    trayectorias[[i + 1]] <- listaVectores
  }
  trayectorias2 <- list() # separamos cada trayectoria por ventana
  for (i in 1:length(trayectorias)) {
    auxiliar <- trayectorias[[i]] # Ordenamos por ventana dentro de la trayectoria
    trayectorias[[i]] <- auxiliar[order(sapply(auxiliar, `[`, 2))]
    v <- trayectorias[[i]][[1]][2] # menor ventana de la trayectoria
    vf <- trayectorias[[i]][[length(trayectorias[[i]])]][2] # mayor ventana
    trayectoria <- list()
    aux <- 1
    while (v <= vf) {
      ventana <- list()
      aux2 <- 1
      for (j in 1:length(trayectorias[[i]])) {
        if (v == trayectorias[[i]][[j]][2]) {
          ventana[[aux2]] <- trayectorias[[i]][[j]]
          aux2 <- aux2 + 1
        }
      }
      if (length(ventana) != 0) {
        trayectoria[[aux]] <- ventana
        aux <- aux + 1
      }
      v <- v + 1
    }
    trayectorias2[[i]] <- trayectoria
  }
  # Transformamos a vector la información de cada ventana
  for (i in 1:length(trayectorias2)) {
    for (j in 1:length(trayectorias2[[i]])) {
      ids <- c()
      names <- c()
      for (k in 1:length(trayectorias2[[i]][[j]])) { # Iteramos cada ventana
        id <- trayectorias2[[i]][[j]][[k]][1] # Extraemos el ID del vértice
        ids <- c(ids, id)
        nombre <- nombres[[id]] # Extraemos el nombre del vértice
        names <- c(names, nombre)
      }
      ven <- trayectorias2[[i]][[j]][[1]][2]
      fec <- fechas[[ven]]
      lista <- list(fec, ven, ids, names)
      trayectorias2[[i]][[j]] <- lista
    }
  }
  return(trayectorias2)
}

# ---- Uso ----
# El ID más alto de comunidad es el valor que `infoAluvial` imprimió en consola.
numeroComunidades <- max(sapply(infoNODOS, function(v) v[3]))
informacionTrayectorias <- comunidadTrayectoria(infoNODOS, numeroComunidades,
                                                dfnombres[[2]], dffechas[[2]])
saveRDS(informacionTrayectorias, file = archivo_salida)
message("Comunidades temporales: ", numeroComunidades,
        " | guardado en ", archivo_salida)
