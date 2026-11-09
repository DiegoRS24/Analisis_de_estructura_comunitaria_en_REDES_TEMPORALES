# =============================================================================
# 05_supervivencia.R
#
# Descripción:
#   Análisis de supervivencia de las comunidades temporales. Cada comunidad
#   temporal (cada "vida", si una comunidad muere y resurge) es un individuo:
#   su tiempo es el número de ventanas que vive y el evento es su muerte. Si
#   sigue viva en la última ventana, el dato es censurado.
#
#   Realiza cinco tareas, después de preparar los datos:
#
#     1. Supervivencia total: curva de Kaplan-Meier y mediana de vida.
#     2. Supervivencia por periodo de nacimiento: curvas de Kaplan-Meier,
#        prueba log-rank y comparaciones por pares con ajuste de Bonferroni.
#     3. Supervivencia por tamaño inicial y por tipo de nacimiento: curvas de
#        Kaplan-Meier y prueba log-rank.
#     4. Modelo de riesgos proporcionales de Cox con periodo, tamaño inicial
#        (en logaritmo) y tipo de nacimiento, y prueba de proporcionalidad.
#     5. Genéricos y comunidades más influyentes: vértices que nunca fueron
#        inactivos y comunidades temporales ordenadas por importancia.
#
# Entradas:
#   - eventos.rds                 : salida de 04_muertes_y_nacimientos.R
#   - informacionTrayectorias.rds : salida de 02_seguimiento_completo.R
#   - NodosNames.csv              : salida de 01_optimizar_particiones.R
#   - Fechas.csv                  : salida de 01_optimizar_particiones.R
#                                   (solo para contar las ventanas)
#
# Salidas:
#   - Resultados impresos en consola (curvas, pruebas y modelo).
#   - figuras/*.png (si `guardar_figuras` es TRUE)
#
# Nota:
#   Las comunidades que ya existían en la ventana 1 (tipo de nacimiento
#   "Ventana 1") nacieron antes del periodo observado, por lo que su tiempo de
#   vida está truncado por la izquierda. Con `excluir_ventana1 <- TRUE` se
#   excluyen del análisis.
#
# Uso:
#   1. Ajustar la sección "Configuración".
#   2. Ejecutar:  source("05_supervivencia.R")
#
# Autor: Diego Antonio Rodríguez Sánchez
# Licencia: MIT (ver LICENSE)
# =============================================================================

# ------------------------------- Configuración -------------------------------
library(survival)
library(survminer)
library(ggplot2)

archivo_eventos      <- "eventos.rds"                  # Salida de 04
archivo_trayectorias <- "informacionTrayectorias.rds"  # Salida de 02
archivo_nodos        <- "NodosNames.csv"               # Salida de 01
archivo_fechas       <- "Fechas.csv"                   # Salida de 01

# Periodos de nacimiento
cortes_periodos    <- NULL  # Ventana donde empieza cada periodo posterior al primero.
                            # Ejemplo: c(100, 200) define tres periodos: ventanas 1-99,
                            # 100-199 y 200-última. NULL: un solo periodo.
etiquetas_periodos <- NULL  # Nombres de los periodos. NULL: "P1", "P2", ...

# Tamaño inicial
cortes_tamano <- NULL       # Límites de las categorías. NULL: terciles.
                            # Ejemplo: c(0, 6, 12, Inf)

excluir_ventana1 <- FALSE   # TRUE: excluye las comunidades presentes en la ventana 1

guardar_figuras <- TRUE     # TRUE: guarda cada gráfica como PNG
carpeta_figuras <- "figuras"

# ------------------------------- Carga de datos ------------------------------
for (f in c(archivo_eventos, archivo_trayectorias, archivo_nodos, archivo_fechas)) {
  if (!file.exists(f)) stop("No se encontró el archivo de entrada: ", f)
}
comunidades    <- readRDS(archivo_eventos)
trayectorias   <- readRDS(archivo_trayectorias)
dfnames        <- read.csv(archivo_nodos, stringsAsFactors = FALSE)
numeroVentanas <- nrow(read.csv(archivo_fechas, stringsAsFactors = FALSE))

# ================================= Funciones =================================
#' Prepara la información de las comunidades temporales para supervivencia
#'
#' @param comunidades  Lista `eventos.rds`.
#' @param n            Número de ventanas de la red.
#' @param periodos     Vector que indica la partición del tiempo en periodos. Si
#'                     hay un periodo de 1 a t y un segundo de t + k a t2, el
#'                     vector es (1, t, t + k, t2); todos valores numéricos.
#' @param trayectorias Lista `informacionTrayectorias.rds`.
#' @return Lista de seis vectores: id, tiempo (ventanas que vivió), evento (1 si
#'         hay muerte registrada, 0 si vive en la última ventana), periodo en
#'         que nació, tipo de nacimiento y tamaño inicial.
supervivencia <- function(comunidades, n, periodos, trayectorias) {
  lista <- list() # Lista a devolver
  auxiliar <- length(comunidades) # número de comunidades temporales detectadas
  id <- 1:auxiliar # vector de ids
  tiempo <- c() # Número de ventanas que vivió
  evento <- c() # 1 si hay muerte registrada; 0 si vive en última ventana
  periodo <- c() # Época en la que nació la comunidad
  tipo_nac <- c() # Clasificación del tipo de nacimiento
  tamano_inicial <- c()
  for (i in 1:auxiliar) {
    c <- comunidades[[i]]
    # tiempo i
    t <- c[2] - c[1] + 1
    tiempo <- c(tiempo, t)
    # evento i
    if (c[2] == n) {
      evento <- c(evento, 0)
    } else {
      evento <- c(evento, 1)
    }
    # periodo i
    p <- 0
    aux <- 2
    while (aux <= length(periodos)) {
      if (c[1] >= periodos[aux - 1] && c[1] <= periodos[aux]) {
        p <- aux / 2
        break
      } else {
        aux <- aux + 2
      }
    }
    periodo <- c(periodo, p)
    # tipo_nac i
    t_n <- c[9]
    tipo_nac <- c(tipo_nac, t_n)
    # tamaño_inicial i
    ti <- 0
    aux <- c[3] + 1
    for (j in 1:length(trayectorias[[aux]])) {
      if (trayectorias[[aux]][[j]][[2]] == c[1]) {
        ti <- length(trayectorias[[aux]][[j]][[3]])
      }
    }
    tamano_inicial <- c(tamano_inicial, ti)
  }
  lista <- list(id, tiempo, evento, periodo, tipo_nac, tamano_inicial)
  return(lista)
}

#' Vértices (genéricos) que nunca se unieron a vértices inactivos
#'
#' @param trayectorias Lista `informacionTrayectorias.rds`.
#' @param nombres      Vector con los nombres de los vértices.
#' @return Lista de listas, una por vértice: nombre, TRUE si el vértice fue
#'         inactivo alguna vez (FALSE en otro caso) y, si lo fue, la ventana en
#'         que lo fue por primera vez.
genericosInfluyentes <- function(trayectorias, nombres) {
  comunidad_0 <- trayectorias[[1]] # comunidad 0
  lista <- list() # Aquí guardamos los genéricos y su valor booleano
  for (i in 1:length(nombres)) {
    lista[[i]] <- list(nombres[[i]], FALSE)
  }
  for (i in 1:length(nombres)) {
    for (j in 1:length(comunidad_0)) {
      if (lista[[i]][[2]]) {
        break
      }
      for (k in 1:length(comunidad_0[[j]][[4]])) {
        if (nombres[[i]] == comunidad_0[[j]][[4]][[k]]) {
          aux <- lista[[i]]
          lista[[i]] <- list(aux[[1]], TRUE, j)
          break
        }
      }
    }
  }
  return(lista)
}

#' Enlista las comunidades temporales por "importancia"
#'
#' La importancia es el número de vértices que tiene la comunidad en cada
#' ventana, sumado a lo largo de su vida. La intersección de los vértices por
#' ventana es el "núcleo" de la comunidad (los que siempre estuvieron en ella);
#' la unión son todos los que estuvieron en ella en alguna ventana.
#'
#' @param comunidades  Lista `eventos.rds`.
#' @param trayectorias Lista `informacionTrayectorias.rds`.
#' @return Lista de listas, una por vida de comunidad temporal:
#'         (ID, importancia, ventana de nacimiento, ventana de muerte, unión de
#'         vértices, intersección de vértices, "Original" o "Resurgimiento k").
importancia <- function(comunidades, trayectorias) {
  auxiliar <- trayectorias
  comunidades_final <- list() # Lista a devolver
  vec <- rep(0, length(trayectorias) - 1)
  for (i in 1:length(comunidades)) {
    aux <- comunidades[[i]]
    id <- aux[[3]] # ID
    ven_nac <- aux[[1]] # Ventana nacimiento
    ven_muer <- aux[[2]] # Ventana muerte
    uni <- c()
    inters <- auxiliar[[id + 1]][[1]][[4]]
    importancia <- 0
    for (j in 1:(ven_muer - ven_nac + 1)) {
      com_ven <- auxiliar[[id + 1]][[j]][[4]]
      uni <- union(uni, com_ven)
      inters <- intersect(inters, com_ven)
      importancia <- importancia + length(com_ven)
      auxiliar[[id + 1]][[j]] <- NA
    }
    if (vec[id] == 0) {
      estado <- "Original"
    } else {
      estado <- paste0("Resurgimiento ", vec[id])
    }
    vec[id] <- vec[id] + 1
    auxiliar[[id + 1]] <- auxiliar[[id + 1]][!is.na(auxiliar[[id + 1]])]
    comunidad <- list(id, importancia, ven_nac, ven_muer, uni, inters, estado)
    comunidades_final[[i]] <- comunidad
  }
  return(comunidades_final)
}

#' Muestra una gráfica y, si `guardar_figuras` es TRUE, la guarda como PNG
#' @param grafica Objeto ggplot o ggsurvplot.
#' @param nombre  Nombre del archivo, sin extensión.
muestra_grafica <- function(grafica, nombre, ancho = 10, alto = 6) {
  print(grafica)
  if (guardar_figuras) {
    dir.create(carpeta_figuras, showWarnings = FALSE, recursive = TRUE)
    g <- if (inherits(grafica, "ggsurvplot")) grafica$plot else grafica
    ggsave(file.path(carpeta_figuras, paste0(nombre, ".png")),
           plot = g, width = ancho, height = alto, dpi = 150)
  }
  invisible(grafica)
}

#' Categoriza el tamaño inicial de las comunidades
#' @param x      Vector de tamaños iniciales.
#' @param cortes Límites de las categorías. NULL: terciles.
#' @return Factor con las categorías.
categoriza_tamano <- function(x, cortes = NULL) {
  if (is.null(cortes)) {
    cortes <- unique(quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE))
  }
  cut(x, breaks = cortes, include.lowest = TRUE, dig.lab = 4)
}

# =========================== Preparación de los datos ========================
# Periodos de nacimiento: pares (inicio, fin) de ventanas
if (is.null(cortes_periodos)) {
  periodos <- c(1, numeroVentanas)
} else {
  if (is.unsorted(cortes_periodos, strictly = TRUE) ||
      any(cortes_periodos < 2) || any(cortes_periodos > numeroVentanas)) {
    stop("`cortes_periodos` debe ser creciente y estar entre 2 y ", numeroVentanas, ".")
  }
  periodos <- as.vector(rbind(c(1, cortes_periodos),
                              c(cortes_periodos - 1, numeroVentanas)))
}
n_periodos <- length(periodos) / 2
if (is.null(etiquetas_periodos)) {
  etiquetas_periodos <- paste0("P", seq_len(n_periodos))
}
if (length(etiquetas_periodos) != n_periodos) {
  stop("`etiquetas_periodos` debe tener ", n_periodos, " elementos.")
}

prueba <- supervivencia(comunidades, numeroVentanas, periodos, trayectorias)
superv <- data.frame(
  id             = prueba[[1]],
  tiempo         = prueba[[2]],
  evento         = prueba[[3]],
  periodo        = prueba[[4]],
  tipo_nac       = prueba[[5]],
  tamano_inicial = prueba[[6]]
)
if (any(superv$tamano_inicial <= 0)) {
  stop("Hay comunidades con tamaño inicial 0: revisa ", archivo_trayectorias, ".")
}

# Variables categóricas con etiquetas
superv$periodo  <- factor(superv$periodo, levels = seq_len(n_periodos),
                          labels = etiquetas_periodos)
superv$tipo_nac <- factor(superv$tipo_nac, levels = c(-1, 0, 1, 2, 3, 4),
                          labels = c("Ventana 1", "División c/muerte", "División c/sup",
                                     "Fusión", "Formación difusa", "Teseo"))
if (anyNA(superv$periodo) || anyNA(superv$tipo_nac)) {
  warning("Hay comunidades sin periodo o sin tipo de nacimiento válido; se excluyen.")
  superv <- superv[complete.cases(superv), ]
}
if (excluir_ventana1) {
  superv <- superv[which(superv$tipo_nac != "Ventana 1"), ]
}
superv$periodo  <- droplevels(superv$periodo)
superv$tipo_nac <- droplevels(superv$tipo_nac)
superv$log_tamano  <- log(superv$tamano_inicial)
superv$tamano_cat  <- categoriza_tamano(superv$tamano_inicial, cortes_tamano)

message("Comunidades temporales: ", nrow(superv),
        " | muertes observadas: ", sum(superv$evento == 1),
        " | censuradas: ", sum(superv$evento == 0))

# ========================== 1. Supervivencia total ===========================
km_total <- survfit(Surv(tiempo, evento) ~ 1, data = superv)
print(km_total)
message("Mediana de vida (ventanas): ", unname(summary(km_total)$table["median"]))
muestra_grafica(
  ggsurvplot(km_total, data = superv,
             conf.int = TRUE, surv.median.line = "hv",
             xlab = "Ventanas de vida",
             ylab = "Probabilidad de supervivencia S(t)",
             title = "Supervivencia general de comunidades",
             palette = "steelblue", ggtheme = theme_bw()),
  "km_total")

# ======================= 2. Supervivencia por periodo ========================
if (nlevels(superv$periodo) >= 2) {
  print(table(superv$periodo))
  km_periodo <- survfit(Surv(tiempo, evento) ~ periodo, data = superv)
  print(km_periodo)
  # Prueba log-rank global y comparaciones por pares (ajuste de Bonferroni)
  print(survdiff(Surv(tiempo, evento) ~ periodo, data = superv))
  print(pairwise_survdiff(Surv(tiempo, evento) ~ periodo, data = superv,
                          p.adjust.method = "bonferroni"))
  muestra_grafica(
    ggsurvplot(km_periodo, data = superv, pval = TRUE,
               xlab = "Ventanas de vida",
               ylab = "Probabilidad de supervivencia S(t)",
               title = "Supervivencia por periodo de nacimiento",
               legend.title = "Periodo",
               legend.labs = levels(superv$periodo),
               palette = "jco", ggtheme = theme_bw()),
    "km_periodo")
} else {
  message("Hay un solo periodo: se omite la comparación (ver `cortes_periodos`).")
}

# ============ 3. Supervivencia por tamaño inicial y tipo de nacimiento =======
# Tamaño inicial
if (nlevels(superv$tamano_cat) >= 2) {
  print(table(superv$tamano_cat))
  km_tamano <- survfit(Surv(tiempo, evento) ~ tamano_cat, data = superv)
  print(km_tamano)
  print(survdiff(Surv(tiempo, evento) ~ tamano_cat, data = superv))
  muestra_grafica(
    ggsurvplot(km_tamano, data = superv, pval = TRUE,
               xlab = "Ventanas de vida",
               ylab = "Probabilidad de supervivencia S(t)",
               title = "Supervivencia por tamaño inicial",
               legend.title = "Tamaño inicial",
               legend.labs = levels(superv$tamano_cat),
               palette = "jco", ggtheme = theme_bw()),
    "km_tamano_inicial")
} else {
  message("El tamaño inicial tiene una sola categoría: se omite (ver `cortes_tamano`).")
}

# Tipo de nacimiento
if (nlevels(superv$tipo_nac) >= 2) {
  print(table(superv$tipo_nac))
  km_tipo <- survfit(Surv(tiempo, evento) ~ tipo_nac, data = superv)
  print(km_tipo)
  print(survdiff(Surv(tiempo, evento) ~ tipo_nac, data = superv))
  print(pairwise_survdiff(Surv(tiempo, evento) ~ tipo_nac, data = superv,
                          p.adjust.method = "bonferroni"))
  muestra_grafica(
    ggsurvplot(km_tipo, data = superv, pval = TRUE,
               xlab = "Ventanas de vida",
               ylab = "Probabilidad de supervivencia S(t)",
               title = "Supervivencia por tipo de nacimiento",
               legend.title = "Tipo de nacimiento",
               legend.labs = levels(superv$tipo_nac),
               palette = "jco", ggtheme = theme_bw()),
    "km_tipo_nacimiento")
} else {
  message("Solo hay un tipo de nacimiento: se omite la comparación.")
}

# =========================== 4. Modelo de Cox ================================
# Covariables: periodo y tipo de nacimiento (si tienen al menos dos categorías)
# y el logaritmo del tamaño inicial (el tamaño es muy asimétrico).
terminos <- "log_tamano"
if (nlevels(superv$periodo)  >= 2) terminos <- c("periodo", terminos)
if (nlevels(superv$tipo_nac) >= 2) terminos <- c(terminos, "tipo_nac")
formula_cox <- as.formula(paste("Surv(tiempo, evento) ~", paste(terminos, collapse = " + ")))

# Categoría de referencia del tipo de nacimiento: la más frecuente
superv_cox <- superv
if (nlevels(superv_cox$tipo_nac) >= 2) {
  superv_cox$tipo_nac <- relevel(superv_cox$tipo_nac,
                                 ref = names(which.max(table(superv_cox$tipo_nac))))
}

modelo_cox <- coxph(formula_cox, data = superv_cox)
print(summary(modelo_cox))

# Supuesto de riesgos proporcionales (prueba de residuos de Schoenfeld).
# Si alguna covariable tiene p < 0.05, el supuesto no se sostiene para ella;
# una opción es estratificar por esa variable con strata().
prueba_ph <- cox.zph(modelo_cox)
print(prueba_ph)

# Razones de riesgo (hazard ratios) con sus intervalos de confianza
tryCatch(
  muestra_grafica(ggforest(modelo_cox, data = superv_cox), "cox_razones_de_riesgo"),
  error = function(e) {
    message("No se pudo generar la gráfica de razones de riesgo: ", conditionMessage(e))
  })

# ============ 5. Genéricos y comunidades más influyentes en el tiempo ========
# Genéricos que nunca fueron inactivos
genericos <- genericosInfluyentes(trayectorias, dfnames[[2]])
nunca_inactivos <- character(0)
for (i in seq_along(genericos)) {
  if (length(genericos[[i]]) != 3) {
    nunca_inactivos <- c(nunca_inactivos, genericos[[i]][[1]])
  }
}
message(length(nunca_inactivos), " vértices nunca fueron inactivos:")
print(nunca_inactivos)

# Comunidades temporales ordenadas de mayor a menor importancia
importantes <- importancia(comunidades, trayectorias)
ordenada_por_importancia <- importantes[order(sapply(importantes, `[[`, 2),
                                             decreasing = TRUE)]
resumen_importantes <- data.frame(
  id                 = sapply(ordenada_por_importancia, `[[`, 1),
  importancia        = sapply(ordenada_por_importancia, `[[`, 2),
  ventana_nacimiento = sapply(ordenada_por_importancia, `[[`, 3),
  ventana_muerte     = sapply(ordenada_por_importancia, `[[`, 4),
  estado             = sapply(ordenada_por_importancia, `[[`, 7)
)
print(head(resumen_importantes, 10))
