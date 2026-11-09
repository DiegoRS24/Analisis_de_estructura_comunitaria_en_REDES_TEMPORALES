# Seguimiento de comunidades en redes temporales (R)

Flujo de análisis en R para estudiar cómo evolucionan las comunidades de una red a lo largo del tiempo. A partir de una partición en comunidades por ventana temporal, el código rastrea comunidades temporales, mide la estabilidad global de la estructura, clasifica los nacimientos y las muertes de comunidades y analiza cuánto sobreviven.

## Contexto

Desarrollado durante mi servicio social en el Instituto de Investigaciones Económicas de la UNAM, aplicado a una red de precios. Los datos y los resultados no se incluyen en este repositorio porque están en preparación para su publicación; la metodología es general y funciona con cualquier red de vértices fijos con una partición en comunidades por ventana.

## Qué hace

1. **Emparejamiento entre ventanas:** para cada comunidad de la ventana *t* calcula su mejor emparejamiento por índice de Jaccard con una comunidad de la ventana *t+1*.
2. **Seguimiento de comunidades temporales:** enlaza las comunidades a lo largo del tiempo y genera la información para un diagrama aluvial.
3. **Métricas globales de estabilidad:** resume el comportamiento de toda la red entre ventanas consecutivas.
4. **Eventos de nacimiento y muerte:** clasifica cómo nace y cómo muere cada comunidad temporal y calcula tasas de variación.
5. **Supervivencia:** prepara los datos de vida de las comunidades para el análisis con `survival` y `survminer`.

## Método

- **Emparejamiento y seguimiento:** índice de Jaccard y emparejamiento de mejor match (Jaccard, 1901; Greene et al., 2010). Una comunidad se considera la misma comunidad temporal entre dos ventanas si supera un umbral de Jaccard y, además, conserva una fracción mínima de su primera aparición medida con el índice de inclusión; este segundo criterio evita el problema del barco de Teseo.
- **Estabilidad global:** una métrica que combina, con media armónica, la precisión y el recall del emparejamiento por Jaccard; una versión ajustada que penaliza los periodos con muchos vértices aislados; la proporción de vértices inactivos; y la información mutua normalizada (NMI) entre particiones consecutivas (Danon et al., 2005; implementada en `igraph`).
- **Eventos:** clasificación con el método GED, basado en índices de inclusión (Bródka et al., 2013). Muertes: absorción o fusión, división, fragmentación y Teseo. Nacimientos: división con muerte del origen, división con supervivencia, fusión, formación difusa y Teseo.
- **Supervivencia:** curvas de Kaplan-Meier, pruebas log-rank con ajuste de Bonferroni y modelo de riesgos proporcionales de Cox sobre la vida de las comunidades temporales.

## Estructura del repositorio

```
.
├── README.md
├── docs/                              documentación adicional
└── R/
    ├── 01_optimizar_particiones.R
    ├── 02_seguimiento_completo.R
    ├── 03_metricas_globales.R
    ├── 04_muertes_y_nacimientos.R
    └── 05_supervivencia.R
```

## Requisitos

R y los paquetes `igraph`, `ggplot2`, `tidyr`, `survival` y `survminer`:

```r
install.packages(c("igraph", "ggplot2", "tidyr", "survival", "survminer"))
```

## Datos de entrada

Un archivo `.rds` con una lista de R que contiene un vector de membresía por ventana temporal (como el que devuelve `igraph::membership()`), todos sobre el mismo conjunto de vértices. Por defecto se llama `lista_vectores_de_membresia.rds`; ajusta el nombre en el script 01.

## Cómo ejecutar
La fundamentación matemática completa está en
[`docs/metodologia/metodologia.tex`](docs/metodologia/metodologia.tex).

### Pipeline de ejecución

Ejecuta los scripts **en este orden**. Cada uno documenta sus entradas y
salidas en el encabezado del propio script.

| # | Script | Sección de la metodología | Qué hace |
|---|--------|---------------------------|----------|
| 1 | `01_optimizar_particiones.R` | §2 Rastreo (parte 1) | Extrae fechas y nombres de vértices; convierte `membership` → vectores de IDs; genera `comunidadesporventana.rds`. |
| 2 | `02_seguimiento_completo.R`     | §2 Rastreo (parte 2) | Emparejamiento Jaccard + inclusión con doble umbral; construye las comunidades temporales. |
| 3 | `03_metricas_globales.R`        | §3 Métricas globales | Calcula precisión, recall, métrica armónica, métrica ajustada, NMI y vértices inactivos. |
| 4 | `04_muertes_y_nacimientos.R`    | §4 Eventos + §6.2 Tasas | Clasifica muertes y nacimientos (absorción, división, fragmentación, fusión) y calcula tasas de variación. |
| 5 | `05_supervivencia.R`            | §5 Supervivencia + §6.1 Comunidades importantes | Kaplan–Meier, log-rank, Bonferroni, Cox multivariado; ranking de comunidades por tamaño×longevidad. |

> ⚠️ Los scripts 2–5 dependen de las salidas de los anteriores. Si corres uno
> fuera de orden, revisa el encabezado del script: indica qué archivos `.rds`
> debe encontrar en la carpeta.

El script 04 calcula los eventos una sola vez y los guarda en `eventos.rds`; después puede saltarse esa primera parte y leer el archivo directamente.

### Nota sobre la indexación

La lista `informacionTrayectorias` empieza por la comunidad 0, que reúne a los vértices que no pertenecen a ninguna comunidad no trivial. Para consultar una comunidad temporal por su ID del diagrama aluvial, usa la entrada `i + 1`: `informacionTrayectorias[[i + 1]]`.

## Parámetros

| Parámetro | Dónde | Valor en el análisis original |
|---|---|---|
| Número de ventanas (`numeroVentanas`) | Script 04 | Depende de los datos; ajústalo |
| Umbral de Jaccard del seguimiento | Scripts 02 y 04 | 0.65 |
| Umbrales α y β de GED (muertes y nacimientos) | Script 04 | 0.7 y 0.8 |

## Referencias

- Bródka, P., Saganowski, S., & Kazienko, P. (2013). GED: The method for group evolution discovery in social networks. *Social Network Analysis and Mining*, 3(1), 1–14. https://doi.org/10.1007/s13278-012-0058-8
- Csárdi, G., & Nepusz, T. (2006). The igraph software package for complex network research. *InterJournal, Complex Systems*, 1695(5), 1–9.
- Danon, L., Diaz-Guilera, A., Duch, J., & Arenas, A. (2005). Comparing community structure identification. *Journal of Statistical Mechanics: Theory and Experiment*, 2005(09), P09008.
- Greene, D., Doyle, D., & Cunningham, P. (2010). Tracking the evolution of communities in dynamic social networks. En *2010 International Conference on Advances in Social Networks Analysis and Mining* (pp. 176–183). IEEE.
- Jaccard, P. (1901). Étude comparative de la distribution florale dans une portion des Alpes et du Jura. *Bulletin de la Société Vaudoise des Sciences Naturelles*, 37(142), 547–579.
- Palla, G., Barabási, A.-L., & Vicsek, T. (2007). Quantifying social group evolution. *Nature*, 446(7136), 664–667. https://doi.org/10.1038/nature05670

## Autor

Diego Antonio Rodríguez Sánchez · Licenciatura en Matemáticas, Facultad de Ciencias, UNAM
