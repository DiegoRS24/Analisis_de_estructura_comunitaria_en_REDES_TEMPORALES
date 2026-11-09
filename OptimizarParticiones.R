# setwd("C:/Users/diego/OneDrive/Escritorio/CodigosFinales")
# setwd("C:/Users/IIEc/Documents/Codigos")

miLista  <- readRDS("membresias_desestacionalizada.rds")
miLista <- miLista[[1]] # Depende del formato
#-----------------------------------------------------------
# Este script realiza tres tareas principales apartir de un objeto lista de 
# vectores de membresia, la lista debe contener un vector por cada ventana.
# 1- Extracción de los datos: Nombre de los vértices como archivo csv, fechas
#   correspondientes a cada ventana como archivo csv.
# 2- Optimizar membresia: Genera el archivo comunidadesporventana.rds lo cual 
#   transorma la información de la partición de una ventana como objeto 
#   membership de igraph a vectores con los que  trabajan todas las funciones 
#   del repositorio.
# 3- Best Match Jaccard: Calcula el mejor emparejamiento por jaccard de una 
#   comunidad en una ventana t con una comunidad en la ventana  t+1 para todas
#   las comunidades en la ventana t. Realiza el calculo para cada par de ventanas
#   en la red. Genera en RDS archivo jaccardstotales.rds y recalls.rds
#   utilizados para generar las métricas del análisis global de la red.
# Utiliza el archivo: lista_vectores_de_membresia.rds BDD dada.

#---------------------- 1- Extracción de los datos --------------------------
# Fechas
categorias <- names(miLista) 
df <- data.frame(x = categorias)
lineas <- c(
  '"","x"',
  paste0('"', seq_along(df$x), '",', df$x)
)
writeLines(lineas, "Fechas329.csv")
# Nombre de los vértices
categorias <- names(miLista[[1]]) 
df <- data.frame(x = categorias)
lineas <- c(
  '"","x"',
  paste0('"', seq_along(df$x), '","', df$x, '"')
)
writeLines(lineas, "NodosNames329.csv")

#---------------------- 2- Optimizar membresia --------------------------
# ------------------------------Funciones--------------------------------------
###
# Función que separa en vectores cada comunidad de cada ventana
# el vector consta de los ids de los vértices que pertenecen a la comunidad, según
# el número de vector dentro de la lista generada es el numero de comunidad 
# definido en dentro del objeto membership. 
# Los ids son tomados del archivo NodosNames_fil_0X.csv (base de datos dada).
# La lista devuelta contiene x vectores, con x el número de comunidades 
# dentro de esa ventana. Recibe:
# - obj: tipo de objeto membership
# - df: Archivo csv con los nombres de los vértices
obtenComunidadesPorVentana <- function(obj, df){
  #Lista de listas, cada sublista es una comunidad
  comunidadesVentana <- list()
  # obteniendo el numero total de comunidades por ventana
  x <-  length(unique(obj))
  ##Guardando cada comunidad en una lista en la lista grande
  cont <- 1
  for(j in 1:x){
    comunidad <- c()
    for(i in 1:length(obj)){
      if(obj[i] == j && !is.na(names(obj[i])) ){
        for(k in 1: length(df[[1]])){
          if( names(obj[i]) == df [k, 2]){
            comunidad <- c(comunidad, df[k, 1])
            break
          }
        }
      }
    }
    if(length(comunidad) != 0){
      comunidadesVentana[[cont]] <- comunidad
      cont <- cont+1
    }
  }
  return(comunidadesVentana)
}
# ----------------------------- USO DE FUNCIONES -----------------------------
# Usamos la función para cada ventana de tiempo
df<- read.csv("NodosNames329.csv", stringsAsFactors = FALSE)
comunidadesporventana <- list()
for(i in 1:length(miLista) ){
  comunidadesporventana[[i]] <- obtenComunidadesPorVentana(miLista[[i]], df)
  print(paste0("ventana ", i, " completada"))
}
# Generamos el archivo tipo rds para correrlo una sola vez
saveRDS(comunidadesporventana, file = "comunidadesporventana.rds")

#---------------------- 3- Best Match Jaccard --------------------------
# ------------------------------Funciones--------------------------------------
###
# Función que devuelve una lista con los indices de Jaccard entre cada pareja de
# comunidades de dos ventanas distintas. Recibe:
# - ventanat: Lista de vectores de comunidades (cada comunidad es un vector)
# - ventanatmas1: Lista de vectores de comunidades
jaccardentreVentanas2 <- function(ventanat, ventanatmas1){
  indices <- list()
  y <- 1
  for(i in 1:length(ventanat)){
    for(j in 1:length((ventanatmas1))){
      x <- length(intersect(ventanat[[i]], ventanatmas1[[j]])) / length(union(ventanat[[i]], ventanatmas1[[j]]))
      v <- c(i, j, x)
      indices[[y]] <- v
      y <- y+1
    }
  }
  indices1 <- list()
  z <- 1
  for(i in 1: length(indices)){
    if(indices[[i]][3] != 0 ){
      indices1[[z]] <- indices[[i]]
      z <- z+1
    }
  }
  return (indices1)
}
###
# Función auxiliar que separa las comparaciones de una comunidad.
# Devuelve una lista de listas. Cada sublista contiene  comparaciones de una 
# comunidad contra todas en t+1. Recibe:
# - indices: Lista de índices de Jaccard entre dos ventanas
separaComunidad <- function(indices){
  keys <- sapply(indices, function(v) v[1])    
  split(indices, keys)                          
}
###
# Función auxiliar que identifica el indice más alto de una comunidad en t, con 
# otra comunidad en t+1. Recibe:
# -indicesSeparados: Lista de índices de Jaccard de una comunidad contra todas
bestMatchFiltro <- function(indicesSeparados){
  indicesAltos <- list()
  y <- 1
  for(i in 1:length(indicesSeparados)){
    x <- 0
    for(j in 1: length(indicesSeparados[[i]])){
      if(x < indicesSeparados[[i]][[j]][3]){
        x <- indicesSeparados[[i]][[j]][3]
        indicesAltos[[y]] <- indicesSeparados[[i]][[j]]
      }
    }
    y <- y+1
  }
  return(indicesAltos)
}
###
# Función que realiza el calculo del índice de similitud por Jaccard más alto de
# una comunidad en la ventana t, con alguna en la ventana t+1.
# Hace uso de las funciones anteriores para todas las comparaciones entre 
# ventanas. Recibe:
# - ventanas: Ventanas de la red en formato de optimizar membresia (vectores)
jaccardbestMatchMasivo <- function(ventanas){
  jaccardscompletos <- list()
  for(i in 1:(length(ventanas)-1)){
    aux <- jaccardentreVentanas2(ventanas[[i]], ventanas[[i+1]])
    aux2 <- separaComunidad(aux)
    jaccardscompletos[[i]] <- bestMatchFiltro(aux2)
    print(paste0("Jaccard v_", i, " -> v_", i+1, " completado"))
  }
  return(jaccardscompletos)
}
###
# Función que hace el calculo del recall, medida utilizada en la métrica
# propuesta. Recibe:
# - ventanas: Ventanas de la red en formato de optimizar membresia (vectores)
jaccardbestMatchMasivoRecall <- function(ventanas){
  jaccardscompletos <- list()
  for(i in 1:(length(ventanas)-1)){
    aux <- jaccardentreVentanas2(ventanas[[i+1]], ventanas[[i]])
    aux2 <- separaComunidad(aux)
    jaccardscompletos[[i]] <- bestMatchFiltro(aux2)
    print(paste0("Recall v_", i+1, " -> v_", i, " completado"))
  }
  return(jaccardscompletos)
}
#----------------------- USO DE FUNCIONES  -----------------------------------
# Calculo de Best match Jacard de t -> t + 1
jaccardstotales <- jaccardbestMatchMasivo(comunidadesporventana)
saveRDS(jaccardstotales, file = "jaccardstotales.rds")
# Calculo de Recalls
recalls <- jaccardbestMatchMasivoRecall(comunidadesporventana)
saveRDS(recalls, file = "recalls.rds")