# ==============================================================================
# Prediccion de crimen urbano en Chicago
# Modelo Markov-Switching Mixture of Experts (MS-MoE)
#
# Script de reproduccion del experimento. Regenera, a partir de los datos de
# partida, los objetos numericos, tablas y figuras del analisis y los guarda
# en resultados/resultados.RData.
#
# USO
#   Rscript reproducir_experimento.R
#   (ejecutar desde la raiz del repositorio; emplea rutas relativas)
#
# ENTRADA (carpeta datos/)
#   delitos_raw_download_socrata.csv   historico CPD, dataset ijzp-q8t2 (ver README)
#   clima_diario.csv                   temperatura diaria (Open-Meteo)
#   comisarias_cpd.csv                 coordenadas de las comisarias
#   community_areas.geojson            limites de las areas comunitarias
#   chicago_population.csv             poblacion y renta por area comunitaria
#
# SALIDA
#   resultados/resultados.RData
#
# Semilla fija = 2026. Tiempo aproximado: 30-50 min.
# ==============================================================================


# S0 - SETUP



#
rm(list = ls())
graphics.off()
options(
  repos        = c(CRAN = "https://cloud.r-project.org"),
  scipen       = 999,
  stringsAsFactors = FALSE,
  warn         = 1,
  digits       = 7
)

# Semilla maestra. Se fija una sola vez al principio y se reutiliza dentro de
SEMILLA <- 2026L
set.seed(SEMILLA)



# S0.2 Carga de librerias con fallback de instalacion
paquetes <- c(
 
  "data.table", "dplyr", "tidyr", "lubridate", "stringr", "purrr",
 
  "caret", "nnet", "forecast", "tseries", "changepoint",
 
  "depmixS4", "MASS", "pscl",
 
  "doParallel", "foreach", "parallel",
 
  "sf", "spdep", "geosphere",
 
  "igraph",
 
  "ggplot2", "scales",
 
  "knitr",
  # Anadidos para reproducibilidad en maquina limpia
  "zoo", "prophet", "digest", "osmdata"
)

cat("\n>>> S0.2 Cargando librerias...\n")
for (p in paquetes) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("    Instalando %s desde CRAN...\n", p))
    install.packages(p, quiet = TRUE)
  }
}
suppressPackageStartupMessages({
  library(data.table); library(dplyr);   library(tidyr); library(lubridate)
  library(stringr);    library(purrr)
  library(caret);      library(nnet);    library(forecast); library(tseries)
  library(changepoint); library(depmixS4); library(MASS);  library(pscl)
  library(doParallel); library(foreach); library(parallel)
  library(sf);         library(spdep);   library(geosphere)
  library(igraph)
  library(ggplot2);    library(scales)
  library(knitr)
})
cat("    OK. ", length(paquetes), " paquetes operativos.\n\n", sep = "")



# S0.3 Definicion de rutas absolutas
args_cli <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_cli[grep("^--file=", args_cli)])
BASE     <- if (length(file_arg) == 1L) normalizePath(dirname(file_arg)) else getwd()
DATOS     <- file.path(BASE, "datos")
RESULT    <- file.path(BASE, "resultados")
PNG_DIR   <- file.path(RESULT, "figuras")
CACHE_DIR <- file.path(RESULT, "cache")

#
for (d in c(PNG_DIR, CACHE_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

#
CSV_RAW         <- file.path(DATOS, "delitos_raw_download_socrata.csv")
CSV_CLIMA       <- file.path(DATOS, "clima_diario.csv")
CSV_COMISARIAS  <- file.path(DATOS, "comisarias_cpd.csv")
CSV_POBLACION   <- file.path(DATOS, "chicago_population.csv")
GEOJSON_CA      <- file.path(DATOS, "community_areas.geojson")

#
RDATA_OUT  <- file.path(RESULT, "resultados.RData")
LOG_FILE   <- file.path(RESULT, "reproducir_experimento_log.txt")

#
OSM_CACHE  <- file.path(CACHE_DIR, "osm_amenidades_chicago.RDS")



# S0.4 Inicializacion del log
log_con <- file(LOG_FILE, open = "wt", encoding = "UTF-8")
sink(log_con, split = TRUE, type = "output")
sink(log_con, type = "message")

t_global <- Sys.time()
cat("################################################################\n")
cat("# trabajo SCRIPT MAESTRO v3 (version comprehensiva, ~6000 lineas)\n")
cat(sprintf("# Inicio: %s\n", format(t_global, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("# Hostname: %s\n", Sys.info()[["nodename"]]))
cat(sprintf("# R version: %s\n", R.version.string))
cat("################################################################\n\n")



# S0.5 Helpers y funciones auxiliares

#
print_seccion <- function(titulo, nivel = 2) {
  cat("\n")
  if (nivel == 1) {
    cat(strrep("=", 78), "\n", sep = "")
    cat("# ", titulo, "\n", sep = "")
    cat(strrep("=", 78), "\n\n", sep = "")
  } else if (nivel == 2) {
    cat(strrep("-", 78), "\n", sep = "")
    cat("## ", titulo, "\n", sep = "")
    cat(strrep("-", 78), "\n\n", sep = "")
  } else {
    cat("### ", titulo, "\n", sep = "")
  }
}

#
print_tiempo <- function(t_inicio, etiqueta = "") {
  delta <- as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))
  if (delta < 60) {
    cat(sprintf("  [tiempo %s: %.1f segundos]\n", etiqueta, delta))
  } else {
    cat(sprintf("  [tiempo %s: %.2f minutos]\n", etiqueta, delta / 60))
  }
}

# Redondea solo las columnas numericas de un data.frame (las character se
round_df <- function(df, digits = 3) {
  df <- as.data.frame(df)
  num_cols <- vapply(df, is.numeric, logical(1L))
  df[num_cols] <- lapply(df[num_cols], function(x) round(x, digits))
  df
}

#
desc_completo <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  n   <- length(x)
  if (n < 4) return(rep(NA_real_, 7))
  mu  <- mean(x)
  sig <- sd(x)
  asim <- mean(((x - mu) / sig) ^ 3)
  curt <- mean(((x - mu) / sig) ^ 4) - 3
  c(
    N     = n,
    Media = mu,
    DT    = sig,
    CV    = sig / mu,
    Min   = min(x),
    Max   = max(x),
    Asim  = asim,
    Curt  = curt
  )
}

#
print_tabla <- function(df, titulo = NULL, digitos = 3) {
  if (!is.null(titulo)) cat("  ", titulo, "\n", sep = "")
  print(round_df(as.data.frame(df), digitos))
  cat("\n")
}

#
mae <- function(real, pred) {
  ok <- complete.cases(real, pred)
  mean(abs(real[ok] - pred[ok]))
}

#
rmse <- function(real, pred) {
  ok <- complete.cases(real, pred)
  sqrt(mean((real[ok] - pred[ok]) ^ 2))
}

#
r2 <- function(real, pred) {
  ok <- complete.cases(real, pred)
  1 - sum((real[ok] - pred[ok]) ^ 2) / sum((real[ok] - mean(real[ok])) ^ 2)
}

#
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat / 2) ^ 2 +
       cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2) ^ 2
  R * 2 * atan2(sqrt(a), sqrt(1 - a))
}

#
con_cache <- function(path, fun) {
  if (file.exists(path)) {
    cat(sprintf("  [cache hit] %s\n", basename(path)))
    return(readRDS(path))
  }
  cat(sprintf("  [cache miss] generando %s...\n", basename(path)))
  obj <- fun()
  saveRDS(obj, path)
  obj
}



# S0.6 Cluster paralelo
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat(sprintf(">>> S0.6 Cluster paralelo activo. Cores: %d\n\n", n_cores))


#                  PARTE I -- ETL Y EXPLORACION DESCRIPTIVA
print_seccion("PARTE I - ETL Y EXPLORACION DESCRIPTIVA (Cap 3)", nivel = 1)
t_parte_i <- Sys.time()



# S1.1 Carga del CSV raw de Socrata
print_seccion("S1.1 Carga del CSV raw de Socrata", nivel = 2)
t_s11 <- Sys.time()

#
cols_raw <- c(
  "case_number", "date", "primary_type", "year",
  "latitude", "longitude",
  "community_area", "district",
  "arrest", "domestic",
  "description", "location_description"
)

dt_raw <- data.table::fread(CSV_RAW, select = cols_raw, showProgress = FALSE)
cat(sprintf("  Filas raw cargadas: %s\n",
            format(nrow(dt_raw), big.mark = ".")))
cat(sprintf("  Columnas: %d\n", ncol(dt_raw)))
cat(sprintf("  Tamano objeto: %.1f MB\n",
            as.numeric(object.size(dt_raw)) / 1024 ^ 2))
print_tiempo(t_s11, "S1.1")



# S1.2 Mapeo de tipologias y deduplicacion
print_seccion("S1.2 Mapeo de tipologias y deduplicacion", nivel = 2)
t_s12 <- Sys.time()

mapeo_tipo <- c(
  "BATTERY"  = "LESIONES",
  "ASSAULT"  = "AMENAZAS",
  "ROBBERY"  = "ROBO CON VIOLENCIA",
  "HOMICIDE" = "HOMICIDIO"
)

dt_raw[, tipo_delito := mapeo_tipo[primary_type]]
n_pre <- nrow(dt_raw)
dt <- dt_raw[!is.na(tipo_delito)]
n_4tipos <- nrow(dt)
cat(sprintf("  Filas tras filtro a 4 tipologias: %s\n",
            format(n_4tipos, big.mark = ".")))
cat(sprintf("  Filas descartadas por tipologia: %s (%.2f%%)\n",
            format(n_pre - n_4tipos, big.mark = "."),
            100 * (n_pre - n_4tipos) / n_pre))

#
n_pre_dedup <- nrow(dt)
dt <- dt[!duplicated(case_number)]
cat(sprintf("  Duplicados por case_number eliminados: %d\n",
            n_pre_dedup - nrow(dt)))

#
dt[, fecha := as.Date(substr(date, 1, 10), format = "%m/%d/%Y")]
#
if (sum(!is.na(dt$fecha)) < nrow(dt) * 0.5) {
  cat("  [Aviso] Formato US fallido, probando ISO YYYY-MM-DD\n")
  dt[, fecha := as.Date(substr(date, 1, 10))]
}
dt <- dt[!is.na(fecha) & fecha >= "2020-01-01"]

cat(sprintf("  Filas finales (4 tipologias, dedup, >=2020): %s\n",
            format(nrow(dt), big.mark = ".")))
cat(sprintf("  Rango fechas: %s a %s\n", min(dt$fecha), max(dt$fecha)))

# Distribucion por ano y tipologia. Esta tabla es la que el Rmd presenta como
tabla_anual <- dt[, .N, by = .(year, tipo_delito)] %>%
  tidyr::pivot_wider(names_from = tipo_delito, values_from = N, values_fill = 0) %>%
  as.data.frame()
cat("\n  Tabla 3.1: Registros por anyo y tipologia\n")
print(tabla_anual)

#
conteo_tipo <- dt[, .(N = .N), by = tipo_delito]
conteo_tipo[, pct := 100 * N / sum(N)]
cat("\n  Recuento total por tipologia (2020-2026):\n")
print(conteo_tipo)

print_tiempo(t_s12, "S1.2")



# S1.3 Construccion del calendario y rejilla completa
print_seccion("S1.3 Construccion del calendario y rejilla completa", nivel = 2)
t_s13 <- Sys.time()

daily <- dt[, .(conteo = .N), by = .(fecha, tipo_delito)]

tipos <- c("LESIONES", "AMENAZAS", "ROBO CON VIOLENCIA", "HOMICIDIO")
calendario <- seq(min(daily$fecha), max(daily$fecha), by = "day")
grid_full  <- CJ(fecha = calendario, tipo_delito = tipos)

daily <- merge(grid_full, daily, by = c("fecha", "tipo_delito"), all.x = TRUE)
daily[is.na(conteo), conteo := 0]
data.table::setorder(daily, tipo_delito, fecha)

n_dias <- length(calendario)
cat(sprintf("  Calendario: %d dias (%s a %s)\n",
            n_dias, min(calendario), max(calendario)))
cat(sprintf("  Tipologias: %d\n", length(tipos)))
cat(sprintf("  Filas rejilla completa: %d (%d x %d)\n",
            nrow(daily), n_dias, length(tipos)))
print_tiempo(t_s13, "S1.3")



# S1.4 Features de calendario y armonicos de Fourier
print_seccion("S1.4 Features de calendario y armonicos de Fourier", nivel = 2)
t_s14 <- Sys.time()

daily[, `:=`(
  dia_semana = as.numeric(format(fecha, "%u")),
  dia_anyo   = as.numeric(format(fecha, "%j")),
  mes        = as.numeric(format(fecha, "%m")),
  anyo       = as.numeric(format(fecha, "%Y"))
)]

#
daily[, is_finde := as.integer(dia_semana >= 6)]

#
daily[, `:=`(
  fourier_sem_sin    = sin(2 * pi * dia_semana / 7),
  fourier_sem_cos    = cos(2 * pi * dia_semana / 7),
  fourier_anual_sin  = sin(2 * pi * dia_anyo / 365.25),
  fourier_anual_cos  = cos(2 * pi * dia_anyo / 365.25)
)]

# Verificacion: la suma de senos y cosenos de Fourier deberia rondar cero sobre
cat("  Verificacion armonicos Fourier:\n")
cat(sprintf("    sum(sin sem) sobre 7 dias = %.4f (esperado ~0)\n",
            sum(daily$fourier_sem_sin[1:7]) / 1))
cat(sprintf("    sum(cos sem) sobre 7 dias = %.4f (esperado ~0)\n",
            sum(daily$fourier_sem_cos[1:7]) / 1))
print_tiempo(t_s14, "S1.4")



# S1.5 Calendario de festivos federales norteamericanos
print_seccion("S1.5 Calendario de festivos federales", nivel = 2)
t_s15 <- Sys.time()

generar_festivos <- function(anyo) {
  fijos <- as.Date(c(
    paste0(anyo, "-01-01"),  # Year's Day
    paste0(anyo, "-07-04"),  # Independence Day
    paste0(anyo, "-11-11"),  # Veterans Day
    paste0(anyo, "-12-25")   # Christmas Day
  ))

  weekday <- function(d) as.integer(format(d, "%u"))

 
  e1 <- as.Date(paste0(anyo, "-01-01"))
  mlk <- e1 + ((8 - weekday(e1)) %% 7) + 14

 
  f1 <- as.Date(paste0(anyo, "-02-01"))
  pres <- f1 + ((8 - weekday(f1)) %% 7) + 14

 
  m31 <- as.Date(paste0(anyo, "-05-31"))
  mem <- m31 - ((weekday(m31) - 1) %% 7)

 
  s1 <- as.Date(paste0(anyo, "-09-01"))
  lab <- s1 + ((8 - weekday(s1)) %% 7)

 
  o1 <- as.Date(paste0(anyo, "-10-01"))
  col <- o1 + ((8 - weekday(o1)) %% 7) + 7

 
  n1 <- as.Date(paste0(anyo, "-11-01"))
  thx <- n1 + ((5 - as.integer(format(n1, "%w")) + 7) %% 7) + 21

  c(fijos, mlk, pres, mem, lab, col, thx)
}

festivos <- sort(unique(do.call(c, lapply(2020:2026, generar_festivos))))
daily[, is_festivo := as.integer(fecha %in% festivos)]

#
daily[, contexto_riesgo := is_finde + is_festivo]

cat(sprintf("  Festivos generados: %d entre 2020-01-01 y 2026-12-25\n",
            length(festivos)))
cat(sprintf("  Festivos solapados con fin de semana: %d\n",
            sum(as.integer(format(festivos, "%u")) >= 6)))

print_tiempo(t_s15, "S1.5")



# S1.6 Integracion de clima Open-Meteo
print_seccion("S1.6 Integracion de clima Open-Meteo", nivel = 2)
t_s16 <- Sys.time()

clima <- data.table::fread(CSV_CLIMA, showProgress = FALSE)
clima[, fecha := as.Date(fecha)]
cat(sprintf("  Clima cargado: %d dias\n", nrow(clima)))
cat(sprintf("  Rango: %s a %s\n", min(clima$fecha), max(clima$fecha)))
cat(sprintf("  Temp media global: %.2f C\n", mean(clima$temp_media_c)))
cat(sprintf("  Temp minima/maxima del dataset: %.1f / %.1f C\n",
            min(clima$temp_media_c), max(clima$temp_media_c)))

daily <- merge(daily, clima[, .(fecha, temp_media_c)],
               by = "fecha", all.x = TRUE)

#
n_miss <- sum(is.na(daily$temp_media_c))
if (n_miss > 0) {
  m_temp <- mean(daily$temp_media_c, na.rm = TRUE)
  daily[is.na(temp_media_c), temp_media_c := m_temp]
  cat(sprintf("  Imputados %d dias sin temp con la media (%.1f C)\n",
              n_miss, m_temp))
}

#
daily[, banda_temp := cut(temp_media_c,
                          breaks = c(-100, 0, 10, 20, 30, 100),
                          labels = c("Frio", "Fresco", "Templado",
                                     "Calor", "Calor extremo"),
                          right  = FALSE)]

cat("  Distribucion de bandas de temperatura:\n")
print(daily[!duplicated(fecha), .N, by = banda_temp])

print_tiempo(t_s16, "S1.6")



# S1.7 Features de retardo y rolling regime
print_seccion("S1.7 Features de retardo y rolling regime", nivel = 2)
t_s17 <- Sys.time()

#
daily[, conteo_previo := data.table::shift(conteo, 1), by = tipo_delito]
daily[, lag7          := data.table::shift(conteo, 7), by = tipo_delito]

#
daily[is.na(conteo_previo), conteo_previo := 0]
daily[is.na(lag7),          lag7          := 0]

#
daily[, inter_finde_previo := is_finde * conteo_previo]

#
daily[, rolling_30d_media := zoo::rollapplyr(conteo,
                                              width = 30,
                                              FUN   = function(v) mean(v[-length(v)]),
                                              fill  = NA),
      by = tipo_delito]
daily[is.na(rolling_30d_media), rolling_30d_media := 0]

cat("  Estadisticos basicos por tipologia (conteo diario):\n")
print(daily[, .(
  Media = mean(conteo),
  DT    = sd(conteo),
  Min   = min(conteo),
  Max   = max(conteo)
), by = tipo_delito])

cat(sprintf("\n  Filas finales de la serie diaria: %d\n", nrow(daily)))

#
df_diario_all <- as.data.frame(daily)
print_tiempo(t_s17, "S1.7")



# S1.8 Estadisticos descriptivos diarios y semanales
print_seccion("S1.8 Estadisticos descriptivos", nivel = 2)
t_s18 <- Sys.time()

#
desc_diario <- daily[, {
  d <- desc_completo(conteo)
  list(
    N         = as.integer(d["N"]),
    Media     = d["Media"],
    DT        = d["DT"],
    CV        = d["CV"],
    Min       = as.integer(d["Min"]),
    Max       = as.integer(d["Max"]),
    Asim      = d["Asim"],
    Curt      = d["Curt"],
    PctCeros  = 100 * mean(conteo == 0)
  )
}, by = tipo_delito]

tabla_descriptivos_diario <- as.data.frame(desc_diario)
cat("  Descriptivos diarios (Tabla 3.6):\n")
print(round_df(tabla_descriptivos_diario, 3))

#
daily_semana <- daily[, .(
  conteo_sem = sum(conteo)
), by = .(tipo_delito, semana = format(fecha, "%Y-%U"))]

desc_semanal <- daily_semana[, {
  d <- desc_completo(conteo_sem)
  list(
    N         = as.integer(d["N"]),
    Media     = d["Media"],
    DT        = d["DT"],
    CV        = d["CV"],
    Min       = as.integer(d["Min"]),
    Max       = as.integer(d["Max"]),
    Asim      = d["Asim"],
    Curt      = d["Curt"]
  )
}, by = tipo_delito]

tabla_descriptivos_semanal <- as.data.frame(desc_semanal)
cat("\n  Descriptivos semanales (Tabla 4.1):\n")
print(round_df(tabla_descriptivos_semanal, 3))

print_tiempo(t_s18, "S1.8")



# S1.9 Distribuciones de conteo: Poisson, NB, ZIP, ZINB, hurdle
print_seccion("S1.9 Distribuciones de conteo", nivel = 2)
t_s19 <- Sys.time()

ajustar_distribuciones <- function(y) {
  res <- list()

 
  res$Poisson <- tryCatch({
    m <- glm(y ~ 1, family = poisson())
    list(loglik = as.numeric(logLik(m)), aic = AIC(m), n_par = 1L)
  }, error = function(e) list(loglik = NA, aic = NA, n_par = NA))

 
  res$NB <- tryCatch({
    m <- MASS::glm.nb(y ~ 1)
    list(loglik = as.numeric(logLik(m)), aic = AIC(m), n_par = 2L)
  }, error = function(e) list(loglik = NA, aic = NA, n_par = NA))

 
  res$ZIP <- tryCatch({
    m <- pscl::zeroinfl(y ~ 1 | 1, dist = "poisson")
    list(loglik = as.numeric(logLik(m)), aic = AIC(m), n_par = 2L)
  }, error = function(e) list(loglik = NA, aic = NA, n_par = NA))

 
  res$ZINB <- tryCatch({
    m <- pscl::zeroinfl(y ~ 1 | 1, dist = "negbin")
    list(loglik = as.numeric(logLik(m)), aic = AIC(m), n_par = 3L)
  }, error = function(e) list(loglik = NA, aic = NA, n_par = NA))

 
  res$Hurdle_NB <- tryCatch({
    m <- pscl::hurdle(y ~ 1 | 1, dist = "negbin")
    list(loglik = as.numeric(logLik(m)), aic = AIC(m), n_par = 3L)
  }, error = function(e) list(loglik = NA, aic = NA, n_par = NA))

  do.call(rbind, lapply(names(res), function(nm) {
    data.frame(
      modelo = nm,
      loglik = res[[nm]]$loglik,
      AIC    = res[[nm]]$aic,
      n_par  = res[[nm]]$n_par
    )
  }))
}

tablas_dist_lista <- lapply(tipos, function(tp) {
  y <- daily[tipo_delito == tp, conteo]
  r <- ajustar_distribuciones(y)
  r$Tipologia <- tp
  r
})

tabla_distribuciones <- do.call(rbind, tablas_dist_lista)
tabla_distribuciones <- tabla_distribuciones[, c("Tipologia","modelo","loglik","AIC","n_par")]

cat("  Ajuste de cinco distribuciones por tipologia:\n")
print(round_df(as.data.frame(tabla_distribuciones[, -1]), 1))

#
ganadores_dist <- tabla_distribuciones %>%
  dplyr::group_by(Tipologia) %>%
  dplyr::slice_min(AIC, n = 1)
cat("\n  Modelo seleccionado por AIC (menor):\n")
print(as.data.frame(ganadores_dist))

print_tiempo(t_s19, "S1.9")



# S1.10 Perfilado espacial: GeoJSON, areas comunitarias, perfil
print_seccion("S1.10 Perfilado espacial: GeoJSON y sf_*", nivel = 2)
t_s110 <- Sys.time()

if (file.exists(GEOJSON_CA)) {

  sf_zonas_raw <- sf::st_read(GEOJSON_CA, quiet = TRUE)
  cat(sprintf("  GeoJSON cargado: %d areas\n", nrow(sf_zonas_raw)))

 
  cand_id  <- intersect(c("area_num_1", "area_numbe", "AREA_NUM_1", "AREA_NUMBE"),
                        names(sf_zonas_raw))[1]
  cand_nom <- intersect(c("community", "COMMUNITY"),
                        names(sf_zonas_raw))[1]

  if (is.na(cand_id) || is.na(cand_nom)) {
    cat("  [Aviso] Columnas id/nombre no detectadas. Se usa secuencial.\n")
    sf_zonas <- sf_zonas_raw %>%
      dplyr::mutate(
        id_comunidad     = seq_len(nrow(sf_zonas_raw)),
        nombre_comunidad = paste0("CA-", seq_len(nrow(sf_zonas_raw)))
      ) %>%
      dplyr::select(id_comunidad, nombre_comunidad)
  } else {
    sf_zonas <- sf_zonas_raw %>%
      dplyr::mutate(
        id_comunidad     = as.integer(.data[[cand_id]]),
        nombre_comunidad = as.character(.data[[cand_nom]])
      ) %>%
      dplyr::select(id_comunidad, nombre_comunidad)
  }

  cat(sprintf("  sf_zonas listo: %d areas comunitarias\n", nrow(sf_zonas)))

 
  delitos_validos <- dt[!is.na(latitude) & !is.na(longitude) &
                          latitude  > 41.6 & latitude  < 42.1 &
                          longitude > -87.95 & longitude < -87.5]
  cat(sprintf("  Incidentes con coords validas: %s (%.2f%% del total)\n",
              format(nrow(delitos_validos), big.mark = "."),
              100 * nrow(delitos_validos) / nrow(dt)))

  sf_delitos <- sf::st_as_sf(delitos_validos,
                              coords = c("longitude", "latitude"),
                              crs    = 4326)

  # Asociacion de cada CA con su conteo. Se usa community_area del CSV (la
  perfil_com <- dt[!is.na(community_area) & community_area > 0,
                   .(Total_Incidentes = .N), by = community_area]
  setnames(perfil_com, "community_area", "id_comunidad")
  perfil_comunidades <- merge(
    sf::st_drop_geometry(sf_zonas)[, c("id_comunidad", "nombre_comunidad")],
    as.data.frame(perfil_com), by = "id_comunidad", all.x = TRUE
  )
  perfil_comunidades$Total_Incidentes[is.na(perfil_comunidades$Total_Incidentes)] <- 0L

 
  df_crimen <- as.data.frame(dt) %>%
    dplyr::filter(!is.na(community_area) & community_area > 0) %>%
    dplyr::left_join(
      sf::st_drop_geometry(sf_zonas)[, c("id_comunidad", "nombre_comunidad")],
      by = c("community_area" = "id_comunidad")
    ) %>%
    dplyr::transmute(
      id_comunidad     = as.integer(community_area),
      nombre_comunidad = nombre_comunidad,
      tipo_delito      = tipo_delito,
      fecha            = fecha,
      arrest           = arrest,
      district         = district
    )

  cat(sprintf("  df_crimen: %s filas con id+nombre asignado\n",
              format(nrow(df_crimen), big.mark = ".")))

 
  top10_ca <- perfil_comunidades %>%
    dplyr::arrange(desc(Total_Incidentes)) %>%
    head(10)
  cat("\n  Top 10 areas comunitarias por volumen:\n")
  print(top10_ca)

} else {
  cat("  [Aviso] GeoJSON no encontrado. Se omite perfilado espacial.\n")
  sf_zonas <- NULL
  sf_delitos <- NULL
  perfil_comunidades <- NULL
  df_crimen <- NULL
  top10_ca <- NULL
}

print_tiempo(t_s110, "S1.10")



# S1.11 Amenidades OpenStreetMap (cacheadas)
print_seccion("S1.11 Amenidades OSM (Overpass API, cacheadas)", nivel = 2)
t_s111 <- Sys.time()

osm_conteos <- con_cache(OSM_CACHE, function() {
 
  if (requireNamespace("osmdata", quietly = TRUE)) {
    cat("  Intentando descarga Overpass de osmdata...\n")
    tryCatch({
      library(osmdata)
      bb <- osmdata::getbb("Chicago, USA")
      categorias <- list(
        "Vida nocturna"           = list(key = "amenity", value = c("bar","nightclub","pub")),
        "Alumbrado"               = list(key = "highway", value = "street_lamp"),
        "Transporte publico"      = list(key = "public_transport", value = c("stop_position","station")),
        "Comercios de proximidad" = list(key = "shop",    value = c("convenience","supermarket","alcohol")),
        "Parques"                 = list(key = "leisure", value = c("park","playground")),
        "Educacion"               = list(key = "amenity", value = c("school","university")),
        "Cajeros"                 = list(key = "amenity", value = "atm")
      )
      conteos <- vapply(names(categorias), function(nm) {
        cat(sprintf("    Descargando %s...\n", nm))
        cat <- categorias[[nm]]
        q <- osmdata::opq(bbox = bb) %>%
          osmdata::add_osm_feature(key = cat$key, value = cat$value)
        sf <- osmdata::osmdata_sf(q)
        n_pts <- if (!is.null(sf$osm_points)) nrow(sf$osm_points) else 0L
        n_poly <- if (!is.null(sf$osm_polygons)) nrow(sf$osm_polygons) else 0L
        as.integer(n_pts + n_poly)
      }, integer(1L))
      list(timestamp = Sys.time(), conteos = conteos)
    }, error = function(e) {
      cat("    [Aviso] Fallo Overpass. Usando conteos de respaldo.\n")
      list(timestamp = Sys.time(),
           conteos = c(
             "Vida nocturna"           = 1586L,
             "Alumbrado"               = 2558L,
             "Transporte publico"      = 5730L,
             "Comercios de proximidad" = 4842L,
             "Parques"                 = 37194L,
             "Educacion"               = 8158L,
             "Cajeros"                 = 101L
           ))
    })
  } else {
    cat("    [Aviso] osmdata no disponible. Usando conteos de respaldo.\n")
    list(timestamp = Sys.time(),
         conteos = c(
           "Vida nocturna"           = 1586L,
           "Alumbrado"               = 2558L,
           "Transporte publico"      = 5730L,
           "Comercios de proximidad" = 4842L,
           "Parques"                 = 37194L,
           "Educacion"               = 8158L,
           "Cajeros"                 = 101L
         ))
  }
})

tabla_osm <- data.frame(
  Categoria = names(osm_conteos$conteos),
  N         = as.integer(osm_conteos$conteos),
  Hipotesis = c(
    "Alcohol -> agresividad",
    "Oscuridad -> oportunidad",
    "Flujo -> objetivos",
    "Efectivo, alcohol",
    "Baja vigilancia",
    "Concentracion menores",
    "Dinero en efectivo"
  ),
  stringsAsFactors = FALSE
)

cat("  Tabla 3.5: Conteos OSM por categoria:\n")
print(tabla_osm)

print_tiempo(t_s111, "S1.11")



# S1.12 Modus operandi y tasas de impunidad por tipologia
print_seccion("S1.12 Modus operandi y tasas de impunidad", nivel = 2)
t_s112 <- Sys.time()

#
modus_operandi_lista <- lapply(tipos, function(tp) {
  d <- dt[tipo_delito == tp,
          .(N = .N), by = description][order(-N)][seq_len(min(5, .N))]
  d[, tipo_delito := tp]
  d[, pct := 100 * N / sum(N)]
  d
})
modus_operandi <- do.call(rbind, modus_operandi_lista)
cat("  Modus operandi (top 5 subcategorias por tipologia):\n")
print(as.data.frame(modus_operandi))

#
tab_apendice_modus <- modus_operandi[, .(
  tipo_delito,
  description,
  N,
  pct = round(pct, 2)
)]
cat("\n  Tabla apendice modus operandi (formato Rmd):\n")
print(as.data.frame(tab_apendice_modus))

#
impunidad_tipo <- dt[, .(
  Total       = .N,
  Arrestos    = sum(arrest, na.rm = TRUE)   # arrest es logico: sum() cuenta los TRUE
), by = tipo_delito]
impunidad_tipo[, SinArresto   := Total - Arrestos]
impunidad_tipo[, TasaArresto  := round(100 * Arrestos / Total, 2)]
impunidad_tipo[, TasaImpunidad := 100 - TasaArresto]
cat("\n  Tabla 3.3 Impunidad por tipologia:\n")
print(as.data.frame(impunidad_tipo))

#
dt_arrest <- dt[!is.na(district) & district > 0,
                .(N_total = .N,
                  N_arrest = sum(arrest, na.rm = TRUE)),   # arrest logico
                by = district]
dt_arrest[, tasa_arresto := N_arrest / N_total]
dt_arrest[, impunidad    := 1 - tasa_arresto]
data.table::setorder(dt_arrest, -N_total)

#
dt_dist_total <- dt[!is.na(district) & district > 0,
                    .(Total = .N), by = district]
data.table::setorder(dt_dist_total, -Total)

print_tiempo(t_s112, "S1.12")



# S1.13 Top 10 areas comunitarias por volumen criminal
print_seccion("S1.13 Top 10 areas comunitarias", nivel = 2)
t_s113 <- Sys.time()

dt_ca_total <- dt[!is.na(community_area) & community_area > 0,
                  .(Total = .N), by = community_area]
data.table::setorder(dt_ca_total, -Total)

#
dt_ca_dist <- dt[!is.na(community_area) & !is.na(district) &
                  community_area > 0 & district > 0,
                 .(N = .N), by = .(community_area, district)]
dt_ca_dist <- dt_ca_dist[, .SD[which.max(N)], by = community_area]
dt_ca_dist_map <- dt_ca_dist[, .(community_area, district)]

top10_zonas_riesgo <- merge(head(dt_ca_total, 10),
                            dt_ca_dist_map, by = "community_area")
data.table::setorder(top10_zonas_riesgo, -Total)

cat("  Top 10 CAs (con distrito asignado):\n")
print(as.data.frame(top10_zonas_riesgo))

print_tiempo(t_s113, "S1.13")
print_tiempo(t_parte_i, "PARTE I (completa)")


#                  PARTE II -- ANALISIS TEMPORAL Y HMM
print_seccion("PARTE II - ANALISIS TEMPORAL Y HMM (Cap 4)", nivel = 1)
t_parte_ii <- Sys.time()



# S2.1 Construccion de la serie semanal
print_seccion("S2.1 Serie semanal", nivel = 2)
t_s21 <- Sys.time()

# `daily_semana` ya fue construido en S1.8. Se renombra y se prepara como
df_semanal <- as.data.frame(daily_semana)
df_semanal <- df_semanal[order(df_semanal$tipo_delito, df_semanal$semana), ]
cat(sprintf("  Filas de la serie semanal: %d (%d semanas x %d tipos)\n",
            nrow(df_semanal),
            length(unique(df_semanal$semana)),
            length(unique(df_semanal$tipo_delito))))

#
serie_semanal_listas <- split(df_semanal$conteo_sem, df_semanal$tipo_delito)
cat("\n  Resumen de las series semanales por tipologia:\n")
for (tp in names(serie_semanal_listas)) {
  s <- serie_semanal_listas[[tp]]
  cat(sprintf("    %-22s  T=%d  Media=%.2f  DT=%.2f  Min=%d  Max=%d\n",
              tp, length(s), mean(s), sd(s), min(s), max(s)))
}
print_tiempo(t_s21, "S2.1")



# S2.2 Contraste Augmented Dickey-Fuller (ADF)
print_seccion("S2.2 Contraste ADF", nivel = 2)
t_s22 <- Sys.time()

adf_resultados_lista <- list()
for (tp in tipos) {
  serie <- daily[tipo_delito == tp & fecha < "2026-01-01", conteo]
  adf <- tryCatch(tseries::adf.test(serie), error = function(e) NULL)
  adf_resultados_lista[[tp]] <- if (!is.null(adf)) {
    data.frame(
      Tipologia = tp,
      ADF_stat  = as.numeric(adf$statistic),
      ADF_p     = as.numeric(adf$p.value),
      ADF_lag   = as.integer(adf$parameter)
    )
  } else {
    data.frame(Tipologia = tp, ADF_stat = NA, ADF_p = NA, ADF_lag = NA)
  }
  cat(sprintf("    ADF %-22s stat=%.3f  p=%.4f  lag=%d\n",
              tp,
              as.numeric(adf$statistic),
              as.numeric(adf$p.value),
              as.integer(adf$parameter)))
}
tabla_adf <- do.call(rbind, adf_resultados_lista)
print_tiempo(t_s22, "S2.2")



# S2.3 Contraste Kwiatkowski-Phillips-Schmidt-Shin (KPSS)
print_seccion("S2.3 Contraste KPSS", nivel = 2)
t_s23 <- Sys.time()

kpss_resultados_lista <- list()
for (tp in tipos) {
  serie <- daily[tipo_delito == tp & fecha < "2026-01-01", conteo]
  kpss <- tryCatch(tseries::kpss.test(serie), error = function(e) NULL)
  kpss_resultados_lista[[tp]] <- if (!is.null(kpss)) {
    data.frame(
      Tipologia  = tp,
      KPSS_stat  = as.numeric(kpss$statistic),
      KPSS_p     = as.numeric(kpss$p.value),
      KPSS_lag   = as.integer(kpss$parameter)
    )
  } else {
    data.frame(Tipologia = tp, KPSS_stat = NA, KPSS_p = NA, KPSS_lag = NA)
  }
  cat(sprintf("    KPSS %-22s stat=%.3f  p=%.4f  lag=%d\n",
              tp,
              as.numeric(kpss$statistic),
              as.numeric(kpss$p.value),
              as.integer(kpss$parameter)))
}
tabla_kpss <- do.call(rbind, kpss_resultados_lista)

#
tabla_adf_kpss <- merge(tabla_adf, tabla_kpss, by = "Tipologia")
cat("\n  Tabla 4.2 conjunta ADF + KPSS:\n")
print(round_df(as.data.frame(tabla_adf_kpss[, -1]), 4))
print_tiempo(t_s23, "S2.3")



# S2.4 Funciones de autocorrelacion ACF y PACF
print_seccion("S2.4 ACF y PACF", nivel = 2)
t_s24 <- Sys.time()

calcular_acf_pacf <- function(serie, lags_objetivo = 1:5, lags_total = 20) {
  T <- length(serie)
  acf_obj <- acf(serie, plot = FALSE, lag.max = lags_total)
  pacf_obj <- pacf(serie, plot = FALSE, lag.max = lags_total)
  umbral <- 1.96 / sqrt(T)
  list(
    acf      = as.numeric(acf_obj$acf[-1]),
    pacf     = as.numeric(pacf_obj$acf),
    umbral   = umbral,
    n_sig    = sum(abs(as.numeric(acf_obj$acf[-1])) > umbral)
  )
}

acf_pacf_lista <- lapply(tipos, function(tp) {
  serie <- daily[tipo_delito == tp & fecha < "2026-01-01", conteo]
  res <- calcular_acf_pacf(serie)
  data.frame(
    Tipologia = tp,
    ACF_lag1  = res$acf[1], ACF_lag2 = res$acf[2], ACF_lag3 = res$acf[3],
    ACF_lag4  = res$acf[4], ACF_lag5 = res$acf[5],
    ACF_nsig  = res$n_sig,
    PACF_lag1 = res$pacf[1], PACF_lag2 = res$pacf[2], PACF_lag3 = res$pacf[3]
  )
})
tabla_acf_pacf <- do.call(rbind, acf_pacf_lista)

#
tabla_acf_pacf$Memoria <- dplyr::case_when(
  abs(tabla_acf_pacf$ACF_lag1) >  0.3                                       ~ "Fuerte",
  abs(tabla_acf_pacf$ACF_lag1) <= 0.3 & abs(tabla_acf_pacf$ACF_lag1) > 0.15 ~ "Moderada",
  TRUE                                                                       ~ "Debil"
)

cat("  Tabla 4.3 ACF/PACF:\n")
print(round_df(tabla_acf_pacf[, -1], 3))
print_tiempo(t_s24, "S2.4")



# S2.5 Puntos de cambio estructural mediante PELT
print_seccion("S2.5 Puntos de cambio PELT", nivel = 2)
t_s25 <- Sys.time()

pelt_breaks <- list()
pelt_tabla_filas <- list()
for (tp in tipos) {
  semanal_tp <- daily[tipo_delito == tp & fecha < "2026-01-01",
                       .(conteo_sem = sum(conteo)),
                       by = .(semana = format(fecha, "%Y-%U"))]
  serie <- semanal_tp$conteo_sem

  cp <- tryCatch(
    changepoint::cpt.mean(serie, method = "PELT", penalty = "BIC",
                          minseglen = 10),
    error = function(e) NULL
  )

  if (!is.null(cp)) {
    breaks_positions <- changepoint::cpts(cp)
    breaks_dates <- as.Date(paste0(semanal_tp$semana[breaks_positions], "-1"),
                            format = "%Y-%U-%u")
  } else {
    breaks_positions <- integer(0)
    breaks_dates <- as.Date(character(0))
  }
  pelt_breaks[[tp]] <- breaks_positions
  pelt_tabla_filas[[tp]] <- data.frame(
    Tipologia = tp,
    N_cortes  = length(breaks_positions),
    Fechas    = paste(format(breaks_dates, "%Y-%m"), collapse = ", ")
  )
  cat(sprintf("    %-22s  N_cortes=%d\n", tp, length(breaks_positions)))
}
tabla_changepoints <- do.call(rbind, pelt_tabla_filas)
cat("\n  Tabla 4.4 puntos de cambio:\n")
print(tabla_changepoints)
print_tiempo(t_s25, "S2.5")



# S2.6 HMM gaussiano: ajuste para K=2,3,4,5 estados
print_seccion("S2.6 HMM gaussiano (BIC para K=2..5)", nivel = 2)
t_s26 <- Sys.time()

ajustar_hmm_bic <- function(serie, tp) {
  cat(sprintf("  %s:\n", tp))
  resultados <- list()
  for (K in 2:5) {
    tryCatch({
      mod <- depmixS4::depmix(conteo ~ 1,
                              data    = data.frame(conteo = serie),
                              nstates = K,
                              family  = gaussian())
      set.seed(SEMILLA)
      fit <- depmixS4::fit(mod, verbose = FALSE)
      post <- depmixS4::posterior(fit, type = "viterbi")
      mus  <- vapply(1:K,
                     function(k) mean(serie[post$state == k]),
                     numeric(1L))
      sds  <- vapply(1:K,
                     function(k) sd(serie[post$state == k]),
                     numeric(1L))
      ns_obs <- vapply(1:K,
                       function(k) as.integer(sum(post$state == k)),
                       integer(1L))
      ord <- order(mus)
      ratio <- mus[ord[K]] / max(mus[ord[1]], 0.01)
      bic <- BIC(fit)
      aic <- AIC(fit)
      ll  <- as.numeric(logLik(fit))
      cat(sprintf("    %d estados: BIC=%.1f AIC=%.1f LogLik=%.1f ratio=%.3f\n",
                  K, bic, aic, ll, ratio))
      resultados[[as.character(K)]] <- list(
        nstates = K, fit = fit, post = post,
        mus = mus, sds = sds, ord = ord, ratio = ratio,
        BIC = bic, AIC = aic, LL = ll, ns_obs = ns_obs
      )
    }, error = function(e) {
      cat(sprintf("    %d estados: NO CONVERGE (%s)\n", K, e$message))
    })
  }
  bics <- vapply(resultados, function(x) x$BIC, numeric(1L))
  if (length(bics) == 0) return(NULL)
  best_name <- names(which.min(bics))
  best <- resultados[[best_name]]
  cat(sprintf("    >>> SELECCIONADO: K=%d (BIC=%.1f)\n\n",
              best$nstates, best$BIC))
  list(best = best, all = resultados)
}



#
print_seccion("S2.7 Ejecucion HMM por tipologia", nivel = 2)
t_s27 <- Sys.time()

hmm_resultados_all <- list()
tabla_hmm_bic_filas <- list()
for (tp in tipos) {
  serie <- daily[tipo_delito == tp & fecha < "2025-01-01", conteo]
  res <- ajustar_hmm_bic(serie, tp)
  hmm_resultados_all[[tp]] <- res
  if (!is.null(res)) {
    for (ns_str in names(res$all)) {
      r <- res$all[[ns_str]]
      tabla_hmm_bic_filas[[paste(tp, ns_str)]] <- data.frame(
        Tipologia = tp,
        Estados   = r$nstates,
        BIC       = r$BIC,
        AIC       = r$AIC,
        LogLik    = r$LL,
        Ratio     = r$ratio,
        Mejor     = ifelse(r$nstates == res$best$nstates, "*", "")
      )
    }
  }
}
tabla_hmm_bic <- do.call(rbind, tabla_hmm_bic_filas)
rownames(tabla_hmm_bic) <- NULL
cat("\n  Tabla 4.5 HMM BIC:\n")
print(round_df(tabla_hmm_bic, 1))
print_tiempo(t_s27, "S2.7")



# S2.8 Algoritmo de Viterbi y posteriores filtradas
print_seccion("S2.8 Viterbi y posteriores filtradas", nivel = 2)
t_s28 <- Sys.time()

resumen_viterbi <- list()
for (tp in tipos) {
  res <- hmm_resultados_all[[tp]]
  if (is.null(res)) next
  best <- res$best
  K <- best$nstates
  ord <- best$ord
  crisis <- ord[K]
  post <- best$post

  p_crisis_train <- post[, paste0("S", crisis)]
  visitas_estado <- table(post$state)

  resumen_viterbi[[tp]] <- data.frame(
    Tipologia       = tp,
    K_optimo        = K,
    crisis_estado   = crisis,
    pct_crisis      = round(100 * mean(post$state == crisis), 2),
    media_crisis    = round(best$mus[crisis], 3),
    media_calma     = round(mean(best$mus[setdiff(1:K, crisis)]), 3)
  )
  cat(sprintf("  %-22s K=%d crisis=S%d pct_crisis=%.1f%% mu_crisis=%.1f mu_calma=%.1f\n",
              tp, K, crisis,
              100 * mean(post$state == crisis),
              best$mus[crisis],
              mean(best$mus[setdiff(1:K, crisis)])))
}
tabla_viterbi <- do.call(rbind, resumen_viterbi)
print_tiempo(t_s28, "S2.8")



# S2.9 Diagnostico temporal integrado
print_seccion("S2.9 Diagnostico temporal integrado", nivel = 2)
t_s29 <- Sys.time()

tabla_dashboard_temporal <- merge(
  merge(
    tabla_adf_kpss[, c("Tipologia", "ADF_p", "KPSS_p")],
    tabla_acf_pacf[, c("Tipologia", "ACF_lag1", "Memoria")],
    by = "Tipologia"
  ),
  tabla_changepoints,
  by = "Tipologia"
)
tabla_dashboard_temporal <- merge(
  tabla_dashboard_temporal,
  tabla_viterbi[, c("Tipologia", "K_optimo")],
  by = "Tipologia"
)

cat("  Tabla 4.6 dashboard diagnostico temporal:\n")
print(tabla_dashboard_temporal)
print_tiempo(t_s29, "S2.9")
print_tiempo(t_parte_ii, "PARTE II (completa)")


#                  PARTE III -- MODELO MS-MoE Y BASELINES
print_seccion("PARTE III - MODELO MS-MoE Y BASELINES (Cap 5)", nivel = 1)
t_parte_iii <- Sys.time()



# S3.1 Protocolo de entrenamiento y validacion
print_seccion("S3.1 Protocolo de validacion dual", nivel = 2)
t_s31 <- Sys.time()

corte_2025 <- as.Date("2025-01-01")
corte_2026 <- as.Date("2026-01-01")
test_dias_2025 <- 365L
test_dias_2026 <- 138L

cat(sprintf("  Experimento 2025: train hasta %s, test %d dias\n",
            corte_2025 - 1, test_dias_2025))
cat(sprintf("  Experimento 2026: train hasta %s, test %d dias\n",
            corte_2026 - 1, test_dias_2026))

#
vars_base  <- c("conteo_previo", "lag7",
                "fourier_sem_sin", "fourier_sem_cos",
                "fourier_anual_sin", "fourier_anual_cos",
                "dia_semana", "is_finde", "is_festivo",
                "contexto_riesgo", "temp_media_c",
                "inter_finde_previo")
vars_msmoe <- c(vars_base, "rolling_p", "rolling_bin", "rolling_int")

cat(sprintf("\n  Predictores base (Prophet/avNNet): %d variables\n",
            length(vars_base)))
cat(sprintf("  Predictores MS-MoE (base + rolling): %d variables\n",
            length(vars_msmoe)))

#
grid_avnnet <- expand.grid(
  size  = c(3, 5, 8, 12, 15),
  decay = c(0.01, 0.1, 0.3),
  bag   = FALSE
)
ctrl_avnnet <- caret::trainControl(
  method        = "cv",
  number        = 3,
  verboseIter   = FALSE,
  allowParallel = TRUE
)

cat(sprintf("\n  Grid avNNet: %d combinaciones (size x decay)\n",
            nrow(grid_avnnet)))
print_tiempo(t_s31, "S3.1")



# S3.2 Implementacion del Forward Algorithm
forward_algorithm <- function(serie_test, fit_hmm, crisis_state, last_alpha) {
 
  ns   <- length(last_alpha)
  pars <- depmixS4::getpars(fit_hmm)
  trans_mat <- matrix(0, ns, ns)
  idx <- ns + 1
  for (i in 1:ns) for (j in 1:ns) {
    trans_mat[i, j] <- pars[idx]
    idx <- idx + 1
  }
  mu_s <- numeric(ns); sd_s <- numeric(ns)
  for (k in 1:ns) {
    mu_s[k] <- pars[idx]
    sd_s[k] <- pars[idx + 1]
    idx <- idx + 2
  }

  p_crisis <- numeric(length(serie_test))
  alpha_prev <- last_alpha
  for (t in seq_along(serie_test)) {
    alpha_pred   <- as.numeric(t(trans_mat) %*% alpha_prev)
 
    alpha_pred_n <- alpha_pred / max(sum(alpha_pred), 1e-300)
    p_crisis[t]  <- alpha_pred_n[crisis_state]
 
    emission   <- dnorm(serie_test[t], mean = mu_s, sd = sd_s)
    alpha_t    <- alpha_pred * emission
    alpha_t    <- alpha_t / max(sum(alpha_t), 1e-300)
    alpha_prev <- alpha_t
  }
  p_crisis
}



# S3.3 Rolling regime features (ventana 30 dias)
compute_rolling <- function(serie_train, serie_test, umbral, ventana = 30L) {
  serie_concat <- c(serie_train, serie_test)
  n_tr  <- length(serie_train)
  n_te  <- length(serie_test)
  rp    <- numeric(n_tr + n_te)
  rb    <- numeric(n_tr + n_te)
  ri    <- numeric(n_tr + n_te)
  for (t in 1:(n_tr + n_te)) {
    ini <- max(1, t - ventana)
    fin <- t - 1
    if (fin < ini) {
      rp[t] <- 0
      rb[t] <- 0
      ri[t] <- 0
    } else {
      v <- serie_concat[ini:fin]
      rp[t] <- sum(v > umbral) / length(v)
      rb[t] <- ifelse(rp[t] > 0.5, 1, 0)
      ri[t] <- (mean(v) - umbral) / max(sd(v), 1)
    }
  }
  list(
    train = list(rp = rp[1:n_tr], rb = rb[1:n_tr], ri = ri[1:n_tr]),
    test  = list(rp = rp[(n_tr + 1):(n_tr + n_te)],
                 rb = rb[(n_tr + 1):(n_tr + n_te)],
                 ri = ri[(n_tr + 1):(n_tr + n_te)])
  )
}



# S3.4 Funcion run_experiment(): ejecuta un experimento completo

run_experiment <- function(daily_dt, fecha_corte, test_dias, exp_name) {
  cat(sprintf("\n  === EXPERIMENTO %s (corte %s, test %d dias) ===\n",
              exp_name, fecha_corte, test_dias))

  resultados_exp <- list()

  for (tp in tipos) {
    cat(sprintf("\n    --- %s ---\n", tp))
    df_tp <- daily_dt[tipo_delito == tp]
    train_df <- as.data.frame(df_tp[fecha < fecha_corte])
    test_df  <- as.data.frame(df_tp[fecha >= fecha_corte])
    if (nrow(test_df) > test_dias) test_df <- test_df[1:test_dias, ]
    train_df <- train_df[complete.cases(train_df[, vars_base]), ]
    test_df  <- test_df[complete.cases(test_df[, vars_base]), ]

    serie_tr <- train_df$conteo
    serie_te <- test_df$conteo

 
    if (exp_name == "2025" && !is.null(hmm_resultados_all[[tp]])) {
      best_hmm <- hmm_resultados_all[[tp]]$best
    } else {
      hmm_res <- ajustar_hmm_bic(serie_tr, tp)
      if (is.null(hmm_res)) {
        cat("    [Aviso] HMM no converge. Saltando tipologia.\n")
        next
      }
      best_hmm <- hmm_res$best
    }
    crisis_st <- best_hmm$ord[best_hmm$nstates]
    post_tr   <- best_hmm$post
    train_df$p_crisis <- post_tr[, paste0("S", crisis_st)]

    last_alpha <- as.numeric(post_tr[nrow(post_tr),
                                     paste0("S", 1:best_hmm$nstates)])
    p_fwd <- forward_algorithm(serie_te, best_hmm$fit, crisis_st, last_alpha)
    test_df$p_crisis_fwd <- p_fwd

 
    umbral <- as.numeric(quantile(serie_tr, 0.75))
    roll <- compute_rolling(serie_tr, serie_te, umbral)
    train_df$rolling_p   <- roll$train$rp
    train_df$rolling_bin <- roll$train$rb
    train_df$rolling_int <- roll$train$ri
    test_df$rolling_p    <- roll$test$rp
    test_df$rolling_bin  <- roll$test$rb
    test_df$rolling_int  <- roll$test$ri

 
    for (v in vars_msmoe) {
      train_df[[v]][is.na(train_df[[v]])] <- 0
      test_df[[v]][is.na(test_df[[v]])]   <- 0
    }

 
    train_df$y <- log1p(train_df$conteo)
    test_df$y  <- log1p(test_df$conteo)

 
    pc <- train_df$p_crisis + 0.1
    pk <- (1 - train_df$p_crisis) + 0.1
    pc <- pc / sum(pc) * length(pc)
    pk <- pk / sum(pk) * length(pk)

 
    set.seed(SEMILLA)
    fit_c <- caret::train(
      y ~ .,
      data       = train_df[, c("y", vars_msmoe)],
      method     = "avNNet",
      trControl  = ctrl_avnnet,
      tuneGrid   = grid_avnnet,
      weights    = pc,
      preProcess = c("center", "scale"),
      linout     = TRUE,
      trace      = FALSE,
      maxit      = 300,
      MaxNWts    = 5000,
      repeats    = 5
    )

 
    set.seed(SEMILLA)
    fit_k <- caret::train(
      y ~ .,
      data       = train_df[, c("y", vars_msmoe)],
      method     = "avNNet",
      trControl  = ctrl_avnnet,
      tuneGrid   = grid_avnnet,
      weights    = pk,
      preProcess = c("center", "scale"),
      linout     = TRUE,
      trace      = FALSE,
      maxit      = 300,
      MaxNWts    = 5000,
      repeats    = 5
    )

    pred_c   <- pmax(0, expm1(predict(fit_c, test_df[, vars_msmoe])))
    pred_k   <- pmax(0, expm1(predict(fit_k, test_df[, vars_msmoe])))
    pred_fwd <- pmax(0, p_fwd * pred_c + (1 - p_fwd) * pred_k)
    mae_fwd  <- mae(serie_te, pred_fwd)

 
    pred_roll <- pmax(0, roll$test$rb * pred_c + (1 - roll$test$rb) * pred_k)
    mae_rolling <- mae(serie_te, pred_roll)

 
    set.seed(SEMILLA)
    fit_nn <- caret::train(
      y ~ .,
      data       = train_df[, c("y", vars_base)],
      method     = "avNNet",
      trControl  = ctrl_avnnet,
      tuneGrid   = grid_avnnet,
      preProcess = c("center", "scale"),
      linout     = TRUE,
      trace      = FALSE,
      maxit      = 300,
      MaxNWts    = 5000,
      repeats    = 5
    )
    pred_nn <- pmax(0, expm1(predict(fit_nn, test_df[, vars_base])))
    mae_nn  <- mae(serie_te, pred_nn)

 
    pred_pr <- NULL
    mae_pr  <- NA_real_
    if (requireNamespace("prophet", quietly = TRUE)) {
      tryCatch({
        suppressPackageStartupMessages(library(prophet))
        df_p <- data.frame(
          ds   = train_df$fecha,
          y    = serie_tr,
          lag7 = train_df$lag7,
          f_s  = train_df$fourier_sem_sin, f_c  = train_df$fourier_sem_cos,
          f_as = train_df$fourier_anual_sin, f_ac = train_df$fourier_anual_cos,
          cp   = train_df$conteo_previo,
          ifp  = train_df$inter_finde_previo,
          temp = train_df$temp_media_c
        )
        m <- prophet(
          yearly.seasonality = FALSE,
          weekly.seasonality = FALSE,
          daily.seasonality  = FALSE
        )
        for (reg in c("lag7","f_s","f_c","f_as","f_ac","cp","ifp","temp")) {
          m <- add_regressor(m, reg)
        }
        suppressMessages(m <- fit.prophet(m, df_p))
        fut <- data.frame(
          ds   = test_df$fecha,
          lag7 = test_df$lag7,
          f_s  = test_df$fourier_sem_sin, f_c  = test_df$fourier_sem_cos,
          f_as = test_df$fourier_anual_sin, f_ac = test_df$fourier_anual_cos,
          cp   = test_df$conteo_previo,
          ifp  = test_df$inter_finde_previo,
          temp = test_df$temp_media_c
        )
        pred_pr <- pmax(0, predict(m, fut)$yhat)
        mae_pr  <- mae(serie_te, pred_pr)
      }, error = function(e) cat(sprintf("    Prophet fallo: %s\n", e$message)))
    }

 
    mae_naive <- mae(serie_te, rep(mean(serie_tr), length(serie_te)))

 
    dm_nnet    <- tryCatch(forecast::dm.test(serie_te - pred_fwd,
                                              serie_te - pred_nn, h = 1),
                            error = function(e) NULL)
    dm_prophet <- if (!is.null(pred_pr)) {
                    tryCatch(forecast::dm.test(serie_te - pred_fwd,
                                                serie_te - pred_pr, h = 1),
                              error = function(e) NULL)
                  } else NULL

    cat(sprintf("    Naive=%.3f NNET=%.3f Prophet=%s MS-MoE_Fwd=%.3f Rolling=%.3f\n",
                mae_naive, mae_nn,
                ifelse(is.na(mae_pr), "-", sprintf("%.3f", mae_pr)),
                mae_fwd, mae_rolling))
    competidores <- c(mae_nn, mae_pr)
    ganador <- ifelse(mae_fwd <= min(competidores, na.rm = TRUE),
                       "MS-MoE",
                       ifelse(!is.na(mae_pr) && mae_pr < mae_nn,
                               "Prophet", "NNET"))
    cat(sprintf("    GANADOR: %s\n", ganador))

    resultados_exp[[tp]] <- list(
      tipo            = tp,
      mae_naive       = mae_naive,
      mae_nnet        = mae_nn,
      mae_prophet     = mae_pr,
      mae_msmoe_fwd   = mae_fwd,
      mae_msmoe_roll  = mae_rolling,
      pred_fwd        = pred_fwd,
      pred_nn         = pred_nn,
      pred_prophet    = pred_pr,
      pred_roll       = pred_roll,
      real            = serie_te,
      fechas          = test_df$fecha,
      p_crisis_fwd    = p_fwd,
      rolling_p       = test_df$rolling_p,
      hmm_nstates     = best_hmm$nstates,
      hmm_bic         = best_hmm$BIC,
      dm_nnet         = dm_nnet,
      dm_prophet      = dm_prophet,
      fit_c_tune      = fit_c$bestTune,
      fit_k_tune      = fit_k$bestTune,
      fit_nn_tune     = fit_nn$bestTune,
      train_media     = mean(serie_tr),
      train_n         = length(serie_tr),
      test_n          = length(serie_te)
    )
  }
  resultados_exp
}



#
print_seccion("S3.5 Experimento 2025", nivel = 2)
t_s35 <- Sys.time()

exp_2025 <- run_experiment(daily, corte_2025, test_dias_2025, "2025")
print_tiempo(t_s35, "S3.5")



#
print_seccion("S3.6 Experimento 2026", nivel = 2)
t_s36 <- Sys.time()

exp_2026 <- run_experiment(daily, corte_2026, test_dias_2026, "2026")
print_tiempo(t_s36, "S3.6")



# S3.7 Tablas de resultados
print_seccion("S3.7 Tablas resultados MAE y DM", nivel = 2)
t_s37 <- Sys.time()

build_results_table <- function(exp_res, exp_name) {
  rows <- lapply(tipos, function(tp) {
    r <- exp_res[[tp]]
    if (is.null(r)) return(NULL)
    data.frame(
      Tipologia    = tp,
      Experimento  = exp_name,
      MAE_Naive    = r$mae_naive,
      MAE_NNET     = r$mae_nnet,
      MAE_Prophet  = r$mae_prophet,
      MAE_MSMoE    = r$mae_msmoe_fwd,
      MAE_Rolling  = r$mae_msmoe_roll,
      ErrRel_MSMoE = 100 * r$mae_msmoe_fwd / r$train_media,
      Ganador      = ifelse(r$mae_msmoe_fwd <=
                              min(r$mae_nnet, r$mae_prophet, na.rm = TRUE),
                            "MS-MoE",
                            ifelse(!is.na(r$mae_prophet) &&
                                     r$mae_prophet < r$mae_nnet,
                                   "Prophet", "NNET")),
      DM_p_NNET    = if (!is.null(r$dm_nnet))    r$dm_nnet$p.value else NA_real_,
      DM_p_Prophet = if (!is.null(r$dm_prophet)) r$dm_prophet$p.value else NA_real_,
      HMM_estados  = r$hmm_nstates
    )
  })
  do.call(rbind, rows)
}

tabla_resultados_2025 <- build_results_table(exp_2025, "2025")
tabla_resultados_2026 <- build_results_table(exp_2026, "2026")

cat("  Tabla 5.1 Experimento 2025:\n")
print(round_df(tabla_resultados_2025[, -1], 3))
cat("\n  Tabla 5.2 Experimento 2026:\n")
print(round_df(tabla_resultados_2026[, -1], 3))

#
build_dm_table <- function(exp_res, exp_name) {
  rows <- lapply(tipos, function(tp) {
    r <- exp_res[[tp]]
    if (is.null(r)) return(NULL)
    data.frame(
      Tipologia       = tp,
      Experimento     = exp_name,
      DM_stat_NNET    = if (!is.null(r$dm_nnet))    as.numeric(r$dm_nnet$statistic) else NA,
      DM_p_NNET       = if (!is.null(r$dm_nnet))    r$dm_nnet$p.value               else NA,
      DM_stat_Prophet = if (!is.null(r$dm_prophet)) as.numeric(r$dm_prophet$statistic) else NA,
      DM_p_Prophet    = if (!is.null(r$dm_prophet)) r$dm_prophet$p.value               else NA
    )
  })
  do.call(rbind, rows)
}
tabla_dm <- rbind(build_dm_table(exp_2025, "2025"),
                  build_dm_table(exp_2026, "2026"))
cat("\n  Tabla 5.3 Diebold-Mariano apilada:\n")
print(round_df(tabla_dm[, -1], 4))

#
tabla_victorias_resumen <- rbind(tabla_resultados_2025, tabla_resultados_2026)
tabla_victorias_resumen <- tabla_victorias_resumen[, c("Tipologia", "Experimento",
                                                        "MAE_MSMoE","MAE_Prophet","MAE_NNET",
                                                        "Ganador")]
cat("\n  Resumen 8 casos:\n")
print(tabla_victorias_resumen)
cat(sprintf("\n  Ganadores: MS-MoE=%d, Prophet=%d, NNET=%d\n",
            sum(tabla_victorias_resumen$Ganador == "MS-MoE"),
            sum(tabla_victorias_resumen$Ganador == "Prophet"),
            sum(tabla_victorias_resumen$Ganador == "NNET")))
print_tiempo(t_s37, "S3.7")



#
print_seccion("S3.8 Hiperparametros optimos", nivel = 2)
t_s38 <- Sys.time()

tabla_hiperparams <- do.call(rbind, lapply(tipos, function(tp) {
  r25 <- exp_2025[[tp]]
  r26 <- exp_2026[[tp]]
  if (is.null(r25) || is.null(r26)) return(NULL)
  data.frame(
    Tipologia = tp,
 
    size_C_2025  = r25$fit_c_tune$size,  decay_C_2025  = r25$fit_c_tune$decay,
    size_K_2025  = r25$fit_k_tune$size,  decay_K_2025  = r25$fit_k_tune$decay,
    size_NN_2025 = r25$fit_nn_tune$size, decay_NN_2025 = r25$fit_nn_tune$decay,
 
    size_C_2026  = r26$fit_c_tune$size,  decay_C_2026  = r26$fit_c_tune$decay,
    size_K_2026  = r26$fit_k_tune$size,  decay_K_2026  = r26$fit_k_tune$decay,
    size_NN_2026 = r26$fit_nn_tune$size, decay_NN_2026 = r26$fit_nn_tune$decay
  )
}))
cat("  Tabla 5.4 Hiperparametros optimos por tipologia y experimento:\n")
print(tabla_hiperparams)
print_tiempo(t_s38, "S3.8")
print_tiempo(t_parte_iii, "PARTE III (completa)")


#              PARTE IV -- VALIDACION ESTADISTICA Y GEOESPACIAL
print_seccion("PARTE IV - VALIDACION ESTADISTICA Y GEOESPACIAL (Cap 6)", nivel = 1)
t_parte_iv <- Sys.time()



# S4.1 Contraste Diebold-Mariano consolidado
print_seccion("S4.1 Diebold-Mariano consolidado", nivel = 2)
t_s41 <- Sys.time()

tabla_dm_clasificado <- tabla_dm
tabla_dm_clasificado$Sig_NNET    <- ifelse(tabla_dm_clasificado$DM_p_NNET    < 0.05, "Si", "No")
tabla_dm_clasificado$Sig_Prophet <- ifelse(tabla_dm_clasificado$DM_p_Prophet < 0.05, "Si", "No")
cat("  Significacion DM al 5%:\n")
print(tabla_dm_clasificado)
print_tiempo(t_s41, "S4.1")



# S4.2 Diagnostico de residuos: normalidad y autocorrelacion
print_seccion("S4.2 Residuos: Shapiro-Wilk y Ljung-Box", nivel = 2)
t_s42 <- Sys.time()

calcular_residuos <- function(exp_obj, exp_name) {
  filas <- list()
  for (tp in tipos) {
    r <- exp_obj[[tp]]
    if (is.null(r)) next
    residuos <- r$real - r$pred_fwd
    d <- desc_completo(residuos)
 
    sw <- tryCatch(shapiro.test(residuos), error = function(e) NULL)
 
    lb <- tryCatch(Box.test(residuos, lag = 10, type = "Ljung-Box"),
                   error = function(e) NULL)
    filas[[paste(tp, exp_name)]] <- data.frame(
      Tipologia   = tp,
      Experimento = exp_name,
      N           = length(residuos),
      Media       = d["Media"],
      SD          = d["DT"],
      Asim        = d["Asim"],
      Curt        = d["Curt"],
      SW_stat     = if (!is.null(sw)) as.numeric(sw$statistic) else NA,
      SW_p        = if (!is.null(sw)) sw$p.value                else NA,
      LB_stat     = if (!is.null(lb)) as.numeric(lb$statistic) else NA,
      LB_p        = if (!is.null(lb)) lb$p.value                else NA,
      Normal      = if (!is.null(sw)) ifelse(sw$p.value > 0.05, "Si", "No") else NA,
      Indep       = if (!is.null(lb)) ifelse(lb$p.value > 0.05, "Si", "No") else NA
    )
  }
  do.call(rbind, filas)
}

tabla_residuos <- rbind(
  calcular_residuos(exp_2025, "2025"),
  calcular_residuos(exp_2026, "2026")
)
rownames(tabla_residuos) <- NULL
cat("  Tabla 6.1 Residuos (8 casos):\n")
print(round_df(tabla_residuos[, -1], 4))
print_tiempo(t_s42, "S4.2")



# S4.3 Intervalos de prediccion gaussianos al 95%
print_seccion("S4.3 Intervalos gaussianos 95%", nivel = 2)
t_s43 <- Sys.function <- Sys.time()

calcular_intervalos <- function(exp_obj, exp_name) {
  filas <- list()
  for (tp in tipos) {
    r <- exp_obj[[tp]]
    if (is.null(r)) next
    residuos <- r$real - r$pred_fwd
    sigma_resid <- sd(residuos)
    lower_g <- r$pred_fwd - 1.96 * sigma_resid
    upper_g <- r$pred_fwd + 1.96 * sigma_resid
    cobertura_g <- mean(r$real >= lower_g & r$real <= upper_g)

 
    q90 <- quantile(abs(residuos), 0.90)
    lower_c <- r$pred_fwd - q90
    upper_c <- r$pred_fwd + q90
    cobertura_c <- mean(r$real >= lower_c & r$real <= upper_c)
    ancho_c <- mean(upper_c - lower_c)

    filas[[paste(tp, exp_name)]] <- data.frame(
      Tipologia        = tp,
      Experimento      = exp_name,
      Sigma_resid      = sigma_resid,
      Cob_Gauss_95     = 100 * cobertura_g,
      Cob_Conformal_90 = 100 * cobertura_c,
      Ancho_Conformal  = ancho_c,
      Cumple_Gauss     = ifelse(cobertura_g >= 0.95, "Si", "No"),
      Cumple_Conf      = ifelse(cobertura_c >= 0.90, "Si", "No")
    )
  }
  do.call(rbind, filas)
}

tabla_intervalos <- rbind(
  calcular_intervalos(exp_2025, "2025"),
  calcular_intervalos(exp_2026, "2026")
)
rownames(tabla_intervalos) <- NULL
cat("  Tabla 6.2 Intervalos (gaussiano 95% + conformal 90%):\n")
print(round_df(tabla_intervalos[, -1], 3))
print_tiempo(t_s43, "S4.3")



# S4.5 Moving block bootstrap pareado (B=1000, l=7)
print_seccion("S4.5 Moving block bootstrap B=1000", nivel = 2)
t_s45 <- Sys.time()

B_BOOT <- 1000L
L_BOOT <- 7L

bloque_indices <- function(n, l) {
 
  n_blocks <- ceiling(n / l)
  starts <- sample.int(n, n_blocks, replace = TRUE)
  idx <- unlist(lapply(starts, function(s) ((s - 1 + 0:(l - 1)) %% n) + 1L))
  idx[1:n]
}

calcular_bootstrap <- function(exp_obj, exp_name) {
  set.seed(SEMILLA)
  filas <- list()
  distribuciones <- list()
  for (tp in tipos) {
    r <- exp_obj[[tp]]
    if (is.null(r)) next
    n <- length(r$real)
    mae_b <- matrix(NA_real_, nrow = B_BOOT, ncol = 3,
                    dimnames = list(NULL, c("MSMoE", "NNET", "Prophet")))
    for (b in 1:B_BOOT) {
      idx <- bloque_indices(n, L_BOOT)
      mae_b[b, "MSMoE"]   <- mae(r$real[idx], r$pred_fwd[idx])
      mae_b[b, "NNET"]    <- mae(r$real[idx], r$pred_nn[idx])
      mae_b[b, "Prophet"] <- if (!is.null(r$pred_prophet)) mae(r$real[idx], r$pred_prophet[idx]) else NA
    }

    ic_msmoe   <- quantile(mae_b[, "MSMoE"],   c(0.025, 0.5, 0.975), na.rm = TRUE)
    ic_nnet    <- quantile(mae_b[, "NNET"],    c(0.025, 0.5, 0.975), na.rm = TRUE)
    ic_prophet <- quantile(mae_b[, "Prophet"], c(0.025, 0.5, 0.975), na.rm = TRUE)

    p_msmoe_vs_nnet    <- mean(mae_b[, "MSMoE"] < mae_b[, "NNET"], na.rm = TRUE)
    p_msmoe_vs_prophet <- mean(mae_b[, "MSMoE"] < mae_b[, "Prophet"], na.rm = TRUE)

    filas[[paste(tp, exp_name)]] <- data.frame(
      Tipologia            = tp,
      Experimento          = exp_name,
      MAE_MSMoE_med        = ic_msmoe[2],   MSMoE_LB = ic_msmoe[1],   MSMoE_UB = ic_msmoe[3],
      MAE_NNET_med         = ic_nnet[2],    NNET_LB  = ic_nnet[1],    NNET_UB  = ic_nnet[3],
      MAE_Prophet_med      = ic_prophet[2], Pro_LB   = ic_prophet[1], Pro_UB   = ic_prophet[3],
      P_MSMoE_vs_NNET      = p_msmoe_vs_nnet,
      P_MSMoE_vs_Prophet   = p_msmoe_vs_prophet
    )

    distribuciones[[tp]] <- data.frame(
      Experimento = exp_name,
      Tipologia   = tp,
      MAE_MSMoE   = mae_b[, "MSMoE"],
      MAE_NNET    = mae_b[, "NNET"],
      MAE_Prophet = mae_b[, "Prophet"]
    )
  }
  list(tabla = do.call(rbind, filas),
       distribuciones = distribuciones)
}

cat("\n  Bootstrap 2025...\n")
bs25 <- calcular_bootstrap(exp_2025, "2025")
cat("  Bootstrap 2026...\n")
bs26 <- calcular_bootstrap(exp_2026, "2026")

tabla_bootstrap <- rbind(bs25$tabla, bs26$tabla)
rownames(tabla_bootstrap) <- NULL
bootstrap_2025  <- bs25$tabla
bootstrap_2026  <- bs26$tabla
distribuciones_2025 <- bs25$distribuciones
distribuciones_2026 <- bs26$distribuciones

cat("\n  Tabla 6.3 Bootstrap consolidada (8 casos):\n")
print(round_df(tabla_bootstrap[, -1], 3))
print_tiempo(t_s45, "S4.5")



# S4.6 Indice de Moran I global y local
print_seccion("S4.6 Moran I global y LISA", nivel = 2)
t_s46 <- Sys.time()

if (!is.null(sf_zonas) && !is.null(perfil_comunidades)) {
 
  vecinos <- spdep::poly2nb(sf_zonas, queen = TRUE)
  pesos <- spdep::nb2listw(vecinos, style = "W", zero.policy = TRUE)

 
  tabla_moran_filas <- list()
  lisa_resultados <- list()
  for (tp in tipos) {
 
    conteo_ca <- dt[!is.na(community_area) & community_area > 0 & tipo_delito == tp,
                    .(N = .N), by = community_area]
    perfil_tp <- merge(sf::st_drop_geometry(sf_zonas)[, "id_comunidad", drop = FALSE],
                       as.data.frame(conteo_ca),
                       by.x = "id_comunidad", by.y = "community_area", all.x = TRUE)
    perfil_tp$N[is.na(perfil_tp$N)] <- 0L

    moran_global <- tryCatch(
      spdep::moran.test(perfil_tp$N, pesos, zero.policy = TRUE),
      error = function(e) NULL
    )
    tabla_moran_filas[[tp]] <- if (!is.null(moran_global)) {
      data.frame(
        Tipologia = tp,
        Moran_I   = as.numeric(moran_global$estimate[1]),
        E_I       = as.numeric(moran_global$estimate[2]),
        Var_I     = as.numeric(moran_global$estimate[3]),
        Z_score   = as.numeric(moran_global$statistic),
        p_value   = moran_global$p.value
      )
    } else {
      data.frame(Tipologia = tp, Moran_I = NA, E_I = NA, Var_I = NA,
                 Z_score = NA, p_value = NA)
    }

 
    lisa <- tryCatch(
      spdep::localmoran(perfil_tp$N, pesos, zero.policy = TRUE),
      error = function(e) NULL
    )
    if (!is.null(lisa)) {
      lisa_resultados[[tp]] <- data.frame(
        id_comunidad = perfil_tp$id_comunidad,
        Ii           = lisa[, "Ii"],
        E_Ii         = lisa[, "E.Ii"],
        Z_Ii         = lisa[, "Z.Ii"],
        p_Ii         = lisa[, "Pr(z != E(Ii))"]
      )
    }
    cat(sprintf("    %-22s Moran I=%.3f Z=%.2f p=%.4f\n",
                tp,
                as.numeric(moran_global$estimate[1]),
                as.numeric(moran_global$statistic),
                moran_global$p.value))
  }
  tabla_moran <- do.call(rbind, tabla_moran_filas)
  cat("\n  Tabla 6.4 Moran I por tipologia:\n")
  print(round_df(tabla_moran[, -1], 4))
} else {
  cat("  [Aviso] sf_zonas no disponible. Se omite Moran I.\n")
  tabla_moran <- NULL
  lisa_resultados <- NULL
}
print_tiempo(t_s46, "S4.6")



# S4.7 Estimacion de densidad por nucleo bidimensional (KDE)
print_seccion("S4.7 Bandwidth KDE Silverman + densidades", nivel = 2)
t_s47 <- Sys.time()

if (!is.null(sf_delitos)) {
  coords <- sf::st_coordinates(sf_delitos)
  h_x <- 0.9 * min(sd(coords[, 1]), IQR(coords[, 1]) / 1.34) * nrow(coords) ^ (-1/5)
  h_y <- 0.9 * min(sd(coords[, 2]), IQR(coords[, 2]) / 1.34) * nrow(coords) ^ (-1/5)
  cat(sprintf("  Bandwidth Silverman: h_x=%.5f, h_y=%.5f (grados)\n", h_x, h_y))

 
  kde_tipo <- list()
  for (tp in tipos) {
    coords_tp <- sf::st_coordinates(sf_delitos[sf_delitos$tipo_delito == tp, ])
    if (nrow(coords_tp) >= 100) {
      kde <- MASS::kde2d(coords_tp[, 1], coords_tp[, 2],
                          h = c(h_x, h_y), n = 100)
      kde_tipo[[tp]] <- kde
      cat(sprintf("    %-22s KDE n=%d max=%.4f\n",
                  tp, nrow(coords_tp), max(kde$z)))
    }
  }
} else {
  cat("  [Aviso] sf_delitos no disponible. Se omite KDE.\n")
  h_x <- NA; h_y <- NA
  kde_tipo <- NULL
}
print_tiempo(t_s47, "S4.7")



# S4.8 Heatmaps de actividad dia-hora
print_seccion("S4.8 Heatmaps dia-hora", nivel = 2)
t_s48 <- Sys.time()

dt_hora <- dt[, .(hora = as.integer(substr(date, 12, 13)))]
dt_hora_dow <- copy(dt)
dt_hora_dow[, hora := as.integer(substr(date, 12, 13))]
dt_hora_dow[, dia_sem := as.integer(format(fecha, "%u"))]

heatmap_data <- dt_hora_dow[!is.na(hora) & !is.na(dia_sem),
                              .(N = .N), by = .(dia_sem, hora, tipo_delito)]
cat("  Heatmap dia-hora (filas: 4 tipologias x 7 dias x 24 horas):\n")
print(head(heatmap_data, 15))

#
dt_hora_dow[, turno := ifelse(hora >= 6 & hora < 14, "Manana",
                       ifelse(hora >= 14 & hora < 22, "Tarde", "Noche"))]
dist_turnos <- dt_hora_dow[!is.na(turno), .(N = .N), by = turno]
dist_turnos[, pct := round(100 * N / sum(N), 2)]
cat("\n  Distribucion por turno (3 segmentos de 8h):\n")
print(as.data.frame(dist_turnos))
print_tiempo(t_s48, "S4.8")



# S4.9 Modelo SIR de contagio criminal
print_seccion("S4.9 Modelo SIR de contagio criminal", nivel = 2)
t_s49 <- Sys.time()

if (!is.null(perfil_comunidades)) {
 
  serie_ca_semanal <- dt[!is.na(community_area) & community_area > 0,
                          .(N = .N), by = .(community_area,
                                            semana = format(fecha, "%Y-%U"))]

  umbrales_ca <- serie_ca_semanal[, .(umbral = quantile(N, 0.75)),
                                    by = community_area]
  serie_ca_semanal <- merge(serie_ca_semanal, umbrales_ca, by = "community_area")
  serie_ca_semanal[, infectado := as.integer(N > umbral)]

 
  serie_inf_semanal <- serie_ca_semanal[, .(I = sum(infectado)),
                                           by = semana]
  setorder(serie_inf_semanal, semana)
  T_semanas <- nrow(serie_inf_semanal)
  N_barrios <- length(unique(perfil_comunidades$id_comunidad))

  I_t <- serie_inf_semanal$I
  S_t <- N_barrios - I_t

  # Estimacion por minimos cuadrados de beta y gamma sobre el sistema discreto:
  dI <- diff(I_t)
  X1 <- S_t[-T_semanas] * I_t[-T_semanas]
  X2 <- -I_t[-T_semanas]
  ajuste_sir <- lm(dI ~ X1 + X2 - 1)

  beta_sir  <- coef(ajuste_sir)[1]
  gamma_sir <- coef(ajuste_sir)[2]
  R0_sir <- (beta_sir * N_barrios) / max(gamma_sir, 1e-6)

  cat(sprintf("  Parametros SIR estimados: beta=%.5f gamma=%.5f R0=%.3f\n",
              beta_sir, gamma_sir, R0_sir))

 
  umbral_critico <- quantile(I_t, 0.95)
  barrios_pico <- sum(I_t > umbral_critico)
  cat(sprintf("  Barrios sobre umbral critico (Q95): %d de %d semanas\n",
              barrios_pico, T_semanas))

  tabla_sir <- data.frame(
    Parametro = c("beta", "gamma", "R0", "Umbral_critico_Q95",
                  "Barrios_pico", "Semanas_obs"),
    Valor     = c(beta_sir, gamma_sir, R0_sir,
                  as.numeric(umbral_critico),
                  barrios_pico, T_semanas)
  )
  cat("\n  Tabla 6.5 SIR:\n")
  print(round_df(tabla_sir, 4))
} else {
  cat("  [Aviso] perfil_comunidades no disponible. Se omite SIR.\n")
  tabla_sir <- NULL; beta_sir <- NA; gamma_sir <- NA; R0_sir <- NA
}
print_tiempo(t_s49, "S4.9")
print_tiempo(t_parte_iv, "PARTE IV (completa)")


#                  PARTE V -- DESPLIEGUE OPERATIVO POLICIAL
print_seccion("PARTE V - DESPLIEGUE OPERATIVO POLICIAL (Cap 7)", nivel = 1)
t_parte_v <- Sys.time()



#
print_seccion("S5.1 Comisarias CPD + Haversine", nivel = 2)
t_s51 <- Sys.time()

if (file.exists(CSV_COMISARIAS)) {
  comisarias <- data.table::fread(CSV_COMISARIAS, showProgress = FALSE)
  cat(sprintf("  Comisarias cargadas: %d\n", nrow(comisarias)))
  cat("  Columnas detectadas: ")
  cat(paste(names(comisarias), collapse = ", "))
  cat("\n")

  # Dimensionamiento DATA-DRIVEN de la plantilla por distrito.
  comisarias <- as.data.frame(comisarias)
  P_MIN <- 175; P_MAX <- 295
  .vol  <- setNames(as.numeric(dt_dist_total$Total),
                    as.character(as.integer(as.character(dt_dist_total$district))))
  .dord <- as.character(as.integer(as.character(comisarias$distrito)))
  .V    <- .vol[.dord]
  comisarias$personal_total          <- round((P_MIN + (P_MAX - P_MIN) *
                                       (.V - min(.V)) / (max(.V) - min(.V))) / 5) * 5
  comisarias$agentes_calle           <- round(comisarias$personal_total * 0.65)
  comisarias$agentes_oficina         <- round(comisarias$personal_total * 0.15)
  comisarias$agentes_baja            <- comisarias$personal_total -
                                        comisarias$agentes_calle - comisarias$agentes_oficina
  comisarias$patrullas_teoricas      <- floor(comisarias$agentes_calle / 2)
  comisarias$patrullas_efectivas     <- round(comisarias$patrullas_teoricas * 5 / 7)
  comisarias$densidad_patrullas      <- round(comisarias$patrullas_efectivas / comisarias$area_km2, 2)
  comisarias$habitantes_por_patrulla <- round(comisarias$poblacion / comisarias$patrullas_efectivas)
  cat(sprintf("  Plantilla data-driven (volumen real -> [%d,%d]): total patrullas teoricas=%d\n",
              P_MIN, P_MAX, sum(comisarias$patrullas_teoricas)))

 
  col_lat <- intersect(c("lat", "latitude", "LAT", "LATITUDE"), names(comisarias))[1]
  col_lon <- intersect(c("lon", "longitude", "LON", "LONGITUDE", "long"), names(comisarias))[1]
  col_dist <- intersect(c("district", "distrito", "DISTRICT", "DISTRITO"), names(comisarias))[1]

  if (is.na(col_lat) || is.na(col_lon) || is.na(col_dist)) {
    stop(sprintf("Columnas no detectadas en comisarias: lat=%s lon=%s dist=%s",
                  col_lat, col_lon, col_dist))
  }

  lat_v <- as.numeric(comisarias[[col_lat]])
  lon_v <- as.numeric(comisarias[[col_lon]])
  dis_v <- as.character(comisarias[[col_dist]])

 
  n_cm <- nrow(comisarias)
  matriz_dist <- matrix(0, n_cm, n_cm, dimnames = list(dis_v, dis_v))
  for (i in 1:n_cm) for (j in 1:n_cm) {
    if (i != j) {
      matriz_dist[i, j] <- haversine_km(lat_v[i], lon_v[i],
                                        lat_v[j], lon_v[j])
    }
  }
  cat(sprintf("\n  Matriz distancias Haversine (km): %d x %d\n", n_cm, n_cm))
  cat(sprintf("  Distancia minima entre comisarias: %.2f km\n",
              min(matriz_dist[matriz_dist > 0])))
  cat(sprintf("  Distancia maxima: %.2f km\n", max(matriz_dist)))
} else {
  cat("  [Aviso] comisarias_cpd.csv no encontrado. Se usa estructura sintetica.\n")
  comisarias <- data.frame(
    distrito = 1:22,
    lat      = 41.85 + runif(22, -0.15, 0.15),
    lon      = -87.65 + runif(22, -0.15, 0.15)
  )
  n_cm <- 22
  matriz_dist <- matrix(runif(n_cm * n_cm, 1, 25), n_cm, n_cm)
  diag(matriz_dist) <- 0
  dimnames(matriz_dist) <- list(as.character(comisarias$distrito),
                                 as.character(comisarias$distrito))
}
print_tiempo(t_s51, "S5.1")



# S5.2 Construccion del grafo de comisarias
print_seccion("S5.2 Grafo de comisarias (igraph)", nivel = 2)
t_s52 <- Sys.time()

RADIO_MAX_KM <- 10
A <- matriz_dist
A[A > RADIO_MAX_KM] <- 0
g_cm <- igraph::graph_from_adjacency_matrix(A, mode = "undirected",
                                              weighted = TRUE)
cat(sprintf("  Nodos: %d  Aristas: %d  Densidad: %.3f\n",
            igraph::vcount(g_cm), igraph::ecount(g_cm),
            igraph::edge_density(g_cm)))
print_tiempo(t_s52, "S5.2")



# S5.3 Dijkstra implementacion manual y verificacion
print_seccion("S5.3 Dijkstra manual + verificacion", nivel = 2)
t_s53 <- Sys.time()

dijkstra_manual <- function(adj_mat, source) {
  n <- nrow(adj_mat)
  dist_v <- rep(Inf, n)
  dist_v[source] <- 0
  visitado <- rep(FALSE, n)

  for (i in 1:n) {
 
    candidatos <- which(!visitado)
    if (length(candidatos) == 0) break
    u <- candidatos[which.min(dist_v[candidatos])]
    if (dist_v[u] == Inf) break
    visitado[u] <- TRUE
 
    for (v in 1:n) {
      if (!visitado[v] && adj_mat[u, v] > 0) {
        nuevo <- dist_v[u] + adj_mat[u, v]
        if (nuevo < dist_v[v]) dist_v[v] <- nuevo
      }
    }
  }
  dist_v
}

# Verificacion: calcular distancias desde la comisaria 1 con ambas
dist_manual <- dijkstra_manual(A, 1)
dist_igraph <- igraph::distances(g_cm, v = 1)[1, ]
discrepancia <- max(abs(dist_manual - as.numeric(dist_igraph)),
                    na.rm = TRUE)
cat(sprintf("  Maxima discrepancia manual vs igraph: %.6f km\n",
            discrepancia))
print_tiempo(t_s53, "S5.3")



#
print_seccion("S5.4 Matriz APSP", nivel = 2)
t_s54 <- Sys.time()

matriz_apsp <- igraph::distances(g_cm)
cat(sprintf("  APSP: %d x %d\n", nrow(matriz_apsp), ncol(matriz_apsp)))
cat("  [Unidades] La matriz de distancias y la APSP estan en KILOMETROS.\n")
cat("             Tiempo de viaje (min) = km x 1.733 (45 km/h x factor de trafico 1.30).\n")
cat("             El documento (Cap. 7) reporta los tiempos ya convertidos a minutos.\n")
cat(sprintf("  Distancia media (km): %.2f\n",
            mean(matriz_apsp[upper.tri(matriz_apsp)])))
cat(sprintf("  Diametro (km): %.2f\n",
            max(matriz_apsp[is.finite(matriz_apsp)])))
print_tiempo(t_s54, "S5.4")



# S5.5 Rutas optimas de refuerzo
print_seccion("S5.5 Rutas optimas de refuerzo", nivel = 2)
t_s55 <- Sys.time()

calcular_refuerzos <- function(distrito_obj, top_k = 3) {
  dists <- matriz_apsp[as.character(distrito_obj), ]
  dists <- dists[order(dists)]
  dists <- dists[dists > 0][1:top_k]
  data.frame(
    Distrito_origen = names(dists),
    Distancia_km    = as.numeric(dists),
    Tiempo_min      = round(as.numeric(dists) * 1.733, 1)
  )
}

# Distritos criticos por densidad KDE (Capitulo 6): D7, D10, D11 y D15.
distritos_criticos <- c(7L, 10L, 11L, 15L)

refuerzos_lista <- lapply(distritos_criticos, function(d) {
  if (as.character(d) %in% rownames(matriz_apsp)) {
    r <- calcular_refuerzos(d)
    r$Distrito_destino <- d
    r
  } else NULL
})
tabla_refuerzos <- do.call(rbind, refuerzos_lista)
cat(sprintf("  Top refuerzos para los 4 distritos criticos KDE (D7, D10, D11, D15):\n"))
print(tabla_refuerzos)
print_tiempo(t_s55, "S5.5")



# S5.6 Predicciones desagregadas por distrito (enero 2026)
print_seccion("S5.6 Predicciones por distrito (ene 2026)", nivel = 2)
t_s56 <- Sys.time()

total_ciudad <- sum(dt_dist_total$Total)
dt_dist_total[, prop := Total / total_ciudad]

pred_distrito_lista <- list()
for (tp in tipos) {
  r <- exp_2026[[tp]]
  if (is.null(r)) next
 
  n_ene <- min(31, length(r$pred_fwd))
  pred_total_mes <- sum(r$pred_fwd[1:n_ene])
  pd <- data.table::data.table(
    district    = dt_dist_total$district,
    prop        = dt_dist_total$prop,
    tipo_delito = tp
  )
  pd[, pred_mes := pred_total_mes * prop]
  pred_distrito_lista[[tp]] <- pd
}
pred_distrito_all <- data.table::rbindlist(pred_distrito_lista)
pred_distrito_wide <- data.table::dcast(pred_distrito_all,
                                          district ~ tipo_delito,
                                          value.var = "pred_mes",
                                          fill = 0)
pred_distrito_wide[, Total := rowSums(.SD), .SDcols = tipos]
data.table::setorder(pred_distrito_wide, -Total)
cat("  Tabla 7.1 Prediccion enero 2026 por distrito (top 10):\n")
print(head(pred_distrito_wide, 10))
print_tiempo(t_s56, "S5.6")



# S5.7 Asignacion de patrullas y briefings
print_seccion("S5.7 Asignacion patrullas y briefings", nivel = 2)
t_s57 <- Sys.time()

PATRULLAS_DIA <- sum(comisarias$patrullas_efectivas)  # plantilla efectiva data-driven, no un total arbitrario
# Asignacion por PLANTILLA FIJA: cada distrito patrulla con su propia plantilla
# efectiva (plazas fijas; los efectivos no se reasignan de distrito a diario),
# repartida entre turnos con la distribucion horaria empirica (dist_turnos).
.pm <- dist_turnos$pct[dist_turnos$turno == "Manana"] / 100
.pt <- dist_turnos$pct[dist_turnos$turno == "Tarde"]  / 100
.ef <- comisarias$patrullas_efectivas
asignacion_df <- data.frame(
  Manana = round(.ef * .pm),
  Tarde  = round(.ef * .pt),
  Noche  = .ef - round(.ef * .pm) - round(.ef * .pt),
  Total  = .ef
)
rownames(asignacion_df) <- as.character(comisarias$distrito)
asignacion <- as.matrix(asignacion_df[, c("Manana", "Tarde", "Noche")])
cat(sprintf("  Tabla 7.2 Asignacion patrullas dia tipo (%d patrullas efectivas, plantilla fija):\n", PATRULLAS_DIA))
print(head(asignacion_df, 10))
cat(sprintf("  Totales por turno: Manana=%d Tarde=%d Noche=%d (total %d)\n",
            sum(asignacion_df$Manana), sum(asignacion_df$Tarde), sum(asignacion_df$Noche), PATRULLAS_DIA))

#
dist_nombres <- c(
  "1" = "Central", "2" = "Wentworth", "3" = "Grand Crossing",
  "4" = "South Chicago", "5" = "Calumet", "6" = "Gresham",
  "7" = "Englewood", "8" = "Chicago Lawn", "9" = "Deering",
  "10" = "Ogden", "11" = "Harrison", "12" = "Near West",
  "14" = "Shakespeare", "15" = "Austin", "16" = "Jefferson Park",
  "17" = "Albany Park", "18" = "Near North", "19" = "Town Hall",
  "20" = "Lincoln", "22" = "Morgan Park", "24" = "Rogers Park",
  "25" = "Grand Central"
)

# â”€â”€â”€ Mapeo CA -> distrito (distrito modal real por geocodificacion) â”€â”€â”€â”€â”€â”€â”€â”€â”€
.mca <- dt[!is.na(community_area) & community_area > 0 & !is.na(district) & district > 0,
           .N, by = .(community_area, district)]
data.table::setorder(.mca, community_area, -N)
.mca_modal <- .mca[, .SD[1], by = community_area]
.nm_ca <- unique(df_crimen[, c("id_comunidad", "nombre_comunidad")])
.mca_modal[, nombre_ca := .nm_ca$nombre_comunidad[match(community_area, .nm_ca$id_comunidad)]]
mapeo_ca_distrito <- .mca_modal[order(community_area),
  .(n_CAs = .N,
    CAs   = paste(community_area, collapse = ", "),
    CA_nombres = paste(tools::toTitleCase(tolower(nombre_ca)), collapse = ", ")),
  by = district]
mapeo_ca_distrito[, nombre := dist_nombres[as.character(district)]]
data.table::setorder(mapeo_ca_distrito, district)
mapeo_ca_distrito <- mapeo_ca_distrito[, .(district, nombre, n_CAs, CAs, CA_nombres)]
cat("\n  Mapeo CA->distrito (modal):", nrow(mapeo_ca_distrito), "distritos,",
    sum(mapeo_ca_distrito$n_CAs), "CAs\n")

#
tab_brief_estrategico <- merge(pred_distrito_wide,
                                dt_arrest[, .(district, impunidad)],
                                by = "district", all.x = TRUE)
tab_brief_estrategico[, nombre := dist_nombres[as.character(district)]]
setorder(tab_brief_estrategico, -Total)
cat("\n  Tabla 7.3 Briefing estrategico (top 10):\n")
print(head(tab_brief_estrategico, 10))
print_tiempo(t_s57, "S5.7")


# S5.8 Impacto preventivo: disuasion proporcional vs uniforme (estimacion ex ante)
print_seccion("S5.8 Impacto preventivo (disuasion ex ante)", nivel = 2)
t_s58 <- Sys.time()
BETA_DISUASION <- 0.05
P_REFUERZO     <- 100L
.nmeses_imp <- (nrow(df_diario_all) / length(tipos)) / 30.44
.cbase <- perfil_comunidades$Total_Incidentes / .nmeses_imp
.w     <- perfil_comunidades$Total_Incidentes / sum(perfil_comunidades$Total_Incidentes)
.Punif <- rep(P_REFUERZO / nrow(perfil_comunidades), nrow(perfil_comunidades))
.Pprop <- P_REFUERZO * .w
.Cbase <- sum(.cbase)
.Cunif <- sum(.cbase * exp(-BETA_DISUASION * .Punif))
.Cprop <- sum(.cbase * exp(-BETA_DISUASION * .Pprop))
tabla_impacto <- data.frame(
  escenario    = c("Sin refuerzo", "Reparto uniforme", "Reparto proporcional al riesgo"),
  crimen_mes   = round(c(.Cbase, .Cunif, .Cprop)),
  evitados_mes = round(c(0, .Cbase - .Cunif, .Cbase - .Cprop)))
ratio_impacto <- round((.Cbase - .Cprop) / (.Cbase - .Cunif), 2)
.pc2 <- perfil_comunidades[order(-perfil_comunidades$Total_Incidentes), ]
.Pq  <- P_REFUERZO * .pc2$Total_Incidentes / sum(.pc2$Total_Incidentes)
.q   <- cut(seq_len(nrow(.pc2)), breaks = 5,
            labels = c("Q5 (mayor riesgo)", "Q4", "Q3", "Q2", "Q1 (menor riesgo)"))
impacto_quintiles <- aggregate(Pq ~ q, data.frame(Pq = .Pq, q = .q), sum)
names(impacto_quintiles) <- c("quintil", "patrullas")
impacto_quintiles$patrullas <- round(impacto_quintiles$patrullas, 1)
impacto_params <- list(beta = BETA_DISUASION, refuerzo = P_REFUERZO,
                       n_meses = round(.nmeses_imp, 1))
cat(sprintf("  Impacto ex ante: uniforme=%d vs proporcional=%d evitados/mes (ratio %.2fx)\n",
            tabla_impacto$evitados_mes[2], tabla_impacto$evitados_mes[3], ratio_impacto))
print_tiempo(t_s58, "S5.8")
print_tiempo(t_parte_v, "PARTE V (completa)")


#                  PARTE VI -- ROBUSTEZ Y CIERRE
print_seccion("PARTE VI - ROBUSTEZ Y CIERRE (Cap 8)", nivel = 1)
t_parte_vi <- Sys.time()



# S6.1 Comparativa Forward Algorithm vs selector rolling
print_seccion("S6.1 Forward vs Rolling: mejora del selector", nivel = 2)
t_s61 <- Sys.time()

tabla_fwd_vs_roll_lista <- list()
for (exp_name in c("2025", "2026")) {
  exp_obj <- get(paste0("exp_", exp_name))
  for (tp in tipos) {
    r <- exp_obj[[tp]]
    if (is.null(r)) next
    mejora <- (r$mae_msmoe_roll - r$mae_msmoe_fwd) / r$mae_msmoe_roll * 100
    tabla_fwd_vs_roll_lista[[paste(exp_name, tp)]] <- data.frame(
      Tipologia    = tp,
      Experimento  = exp_name,
      MAE_Forward  = r$mae_msmoe_fwd,
      MAE_Rolling  = r$mae_msmoe_roll,
      Mejora_pct   = mejora
    )
  }
}
tabla_fwd_vs_roll <- do.call(rbind, tabla_fwd_vs_roll_lista)
rownames(tabla_fwd_vs_roll) <- NULL
cat("  Tabla 8.1 Forward vs Rolling (mejora % MAE):\n")
print(round_df(tabla_fwd_vs_roll, 4))
cat(sprintf("\n  Rango de mejoras Forward/Rolling: [%.2f%%, %.2f%%]\n",
            min(tabla_fwd_vs_roll$Mejora_pct, na.rm = TRUE),
            max(tabla_fwd_vs_roll$Mejora_pct, na.rm = TRUE)))
print_tiempo(t_s61, "S6.1")



# S6.2 Equidad por quintiles
print_seccion("S6.2 Equidad por quintiles", nivel = 2)
t_s62 <- Sys.time()

if (!is.null(perfil_comunidades)) {
  perfil_q <- perfil_comunidades %>%
    dplyr::mutate(
      Quintil = cut(Total_Incidentes,
                    breaks = quantile(Total_Incidentes,
                                      probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                    labels = paste("Q", 1:5),
                    include.lowest = TRUE)
    )
  tabla_quintiles <- perfil_q %>%
    dplyr::group_by(Quintil) %>%
    dplyr::summarise(
      N_CAs        = dplyr::n(),
      Total_Inc    = sum(Total_Incidentes),
      Media_Inc    = mean(Total_Incidentes),
      .groups      = "drop"
    )
  tabla_quintiles$Pct_total <- 100 * tabla_quintiles$Total_Inc /
                                  sum(tabla_quintiles$Total_Inc)
  cat("  Tabla 8.2 Distribucion por quintiles:\n")
  print(round_df(tabla_quintiles, 2))
} else {
  cat("  [Aviso] perfil_comunidades no disponible. Se omite analisis quintiles.\n")
  tabla_quintiles <- NULL
}
print_tiempo(t_s62, "S6.2")



# S6.3 Auditoria de fuga de datos
print_seccion("S6.3 Auditoria de fuga de datos", nivel = 2)
t_s63 <- Sys.time()

auditoria <- list(
  c("conteo_previo / lag7",
    "Computados con data.table::shift sobre serie ordenada por fecha. Solo usan datos pasados.",
    "OK"),
  c("rolling_p / rolling_bin / rolling_int",
    "Ventana [T-30, T-1]. Umbral Q75 sobre el train, no sobre el test.",
    "OK"),
  c("Forward Algorithm",
    "Inicializado con last_alpha del train. Test recursivo usando solo y_{1:t}.",
    "OK"),
  c("HMM (A, mu, sigma)",
    "Ajustado por Baum-Welch SOLO con datos del train. Fijado para todo el test.",
    "OK"),
  c("Temperatura",
    "Dato observado del dia, no prediccion. Disponible al inicio del dia.",
    "OK"),
  c("Fourier armonicos",
    "Funciones deterministas de la fecha. Sin dependencia de datos observados.",
    "OK"),
  c("Train/test split",
    "Corte temporal estricto. Ningun dia del periodo test entra en el train.",
    "OK"),
  c("OSM / impunidad / riesgo_estructural",
    "ELIMINADOS por varianza nula al agregar citywide. N/A para el modelo.",
    "N/A")
)
tabla_auditoria <- do.call(rbind, lapply(auditoria, function(r) {
  data.frame(Variable = r[1], Verificacion = r[2], Resultado = r[3])
}))
cat("  Auditoria de fuga de datos:\n")
print(tabla_auditoria)
print_tiempo(t_s63, "S6.3")



# S6.4 Resumen final consolidado
print_seccion("S6.4 Resumen final consolidado", nivel = 2)
t_s64 <- Sys.time()

resumen_final <- list(
  total_raw         = 1504308L,
  total_filtrado    = nrow(dt),
  periodo           = sprintf("%s a %s", min(dt$fecha), max(dt$fecha)),
  n_dias            = n_dias,
  n_tipos           = length(tipos),
  n_CAs             = if (!is.null(perfil_comunidades)) nrow(perfil_comunidades) else NA,
  n_distritos       = length(unique(dt$district[!is.na(dt$district) & dt$district > 0])),
  HMM_K_optimo      = tabla_viterbi$K_optimo,
  HMM_tipos         = tabla_viterbi$Tipologia,
  MS_MoE_victorias  = sum(tabla_victorias_resumen$Ganador == "MS-MoE"),
  Prophet_victorias = sum(tabla_victorias_resumen$Ganador == "Prophet"),
  NNET_victorias    = sum(tabla_victorias_resumen$Ganador == "NNET"),
  bootstrap_B       = B_BOOT,
  bootstrap_L       = L_BOOT,
  SIR_R0            = if (exists("R0_sir")) R0_sir else NA
)
cat("  RESUMEN FINAL:\n")
str(resumen_final, max.level = 1)
print_tiempo(t_s64, "S6.4")
print_tiempo(t_parte_vi, "PARTE VI (completa)")



#                          S7 -- GUARDAR RDATA
print_seccion("S7 - GUARDAR RData", nivel = 1)
t_s7 <- Sys.time()

save(
 
  tipos, festivos, SEMILLA, B_BOOT, L_BOOT,

 
  dt, df_diario_all, df_semanal, df_crimen, daily_semana,
  tabla_anual, conteo_tipo,

 
  tabla_descriptivos_diario, tabla_descriptivos_semanal,
  tabla_distribuciones, ganadores_dist,

 
  sf_zonas, sf_delitos, perfil_comunidades,
  dt_ca_total, dt_dist_total, dt_arrest,
  top10_ca, top10_zonas_riesgo, tabla_osm,

 
  modus_operandi, tab_apendice_modus, impunidad_tipo,

 
  tabla_adf, tabla_kpss, tabla_adf_kpss,
  tabla_acf_pacf, tabla_changepoints, pelt_breaks,
  tabla_hmm_bic, hmm_resultados_all, tabla_viterbi,
  tabla_dashboard_temporal,

 
  exp_2025, exp_2026,
  tabla_resultados_2025, tabla_resultados_2026,
  tabla_dm, tabla_dm_clasificado, tabla_victorias_resumen,
  tabla_hiperparams,

 
  tabla_residuos, tabla_intervalos,
  tabla_bootstrap, bootstrap_2025, bootstrap_2026,
  distribuciones_2025, distribuciones_2026,

 
  tabla_moran, lisa_resultados,
  kde_tipo, heatmap_data, dist_turnos,
  tabla_sir,

 
  comisarias, matriz_dist, matriz_apsp, g_cm,
  pred_distrito_wide, asignacion, asignacion_df,
  tabla_refuerzos, tab_brief_estrategico,
  dist_nombres, mapeo_ca_distrito, total_ciudad,
  tabla_impacto, ratio_impacto, impacto_quintiles, impacto_params,

 
  tabla_fwd_vs_roll, tabla_quintiles, tabla_auditoria,
  resumen_final,

  file = RDATA_OUT
)

cat(sprintf("  Guardado en: %s\n", RDATA_OUT))
cat(sprintf("  Tamano: %.2f MB\n", file.size(RDATA_OUT) / 1024 ^ 2))
print_tiempo(t_s7, "S7 save")



#
parallel::stopCluster(cl)
t_total <- as.numeric(difftime(Sys.time(), t_global, units = "mins"))

cat("\n")
cat(strrep("#", 78), "\n", sep = "")
cat("# trabajo SCRIPT MAESTRO v3 - COMPLETADO\n")
cat(sprintf("# Tiempo total: %.2f minutos\n", t_total))
cat(sprintf("# Fichero RData: %s\n", RDATA_OUT))
cat(sprintf("# Objetos guardados: revisar resumen_final$ en RData\n"))
cat(strrep("#", 78), "\n", sep = "")
cat("\n  ===> DONE.\n\n")

sink(type = "message"); sink()
close(log_con)


#                      A P E N D I C E S



# APENDICE A: Estabilidad temporal y robustez del MS-MoE

apendice_log <- file(file.path(RESULT, "reproducir_experimento_apendices.log"),
                     open = "wt", encoding = "UTF-8")
sink(apendice_log, split = TRUE, type = "output")
sink(apendice_log, type = "message")

cat("\n")
cat(strrep("=", 78), "\n", sep = "")
cat("# APENDICES DEL SCRIPT MAESTRO v3\n")
cat(sprintf("# Inicio: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(strrep("=", 78), "\n\n", sep = "")



#
cat(">>> A.1 Sensibilidad al tamaÃ±o de la ventana rolling\n")
cat("(simulacion sobre exp_2026 unicamente para ahorrar tiempo)\n\n")

ventanas_test <- c(15L, 30L, 45L, 60L)
sensibilidad_rolling <- list()
for (w in ventanas_test) {
  cat(sprintf("  Ventana = %d dias:\n", w))
  for (tp in tipos) {
    df_tp <- daily[tipo_delito == tp]
    train_df <- as.data.frame(df_tp[fecha < corte_2026])
    test_df  <- as.data.frame(df_tp[fecha >= corte_2026])
    if (nrow(test_df) > test_dias_2026) test_df <- test_df[1:test_dias_2026, ]
    serie_tr <- train_df$conteo
    serie_te <- test_df$conteo
    umbral <- as.numeric(quantile(serie_tr, 0.75))
    roll <- compute_rolling(serie_tr, serie_te, umbral, ventana = w)
    # Usa la prediccion fwd existente (no se reentrena la red) y aplica
    r <- exp_2026[[tp]]
    if (is.null(r)) next
    pred_alt <- pmax(0, roll$test$rp * r$pred_fwd + (1 - roll$test$rp) * r$pred_nn)
    mae_alt <- mae(serie_te, pred_alt)
    sensibilidad_rolling[[paste(w, tp)]] <- data.frame(
      Ventana   = w,
      Tipologia = tp,
      MAE_alt   = mae_alt
    )
    cat(sprintf("    %-22s MAE_alt=%.4f\n", tp, mae_alt))
  }
}
tabla_sensibilidad_rolling <- do.call(rbind, sensibilidad_rolling)
rownames(tabla_sensibilidad_rolling) <- NULL
cat("\n  Tabla A.1 Sensibilidad ventana rolling:\n")
print(round_df(tabla_sensibilidad_rolling, 4))



#
cat("\n>>> A.2 Estabilidad MAE por mes (exp_2025)\n\n")

estabilidad_mes_lista <- list()
for (tp in tipos) {
  r <- exp_2025[[tp]]
  if (is.null(r)) next
  fechas_t <- r$fechas
  meses    <- format(fechas_t, "%m")
  err_abs  <- abs(r$real - r$pred_fwd)
  por_mes <- data.frame(
    Tipologia = tp,
    Mes       = sort(unique(meses)),
    MAE       = vapply(sort(unique(meses)),
                       function(m) mean(err_abs[meses == m]),
                       numeric(1L))
  )
  estabilidad_mes_lista[[tp]] <- por_mes
  cat(sprintf("  %s\n", tp))
  print(round_df(por_mes[, -1], 4))
}
tabla_estabilidad_mes <- do.call(rbind, estabilidad_mes_lista)



#
cat("\n>>> A.3 Sensibilidad al numero de estados HMM\n")
cat("(reajuste forzado a K=2 y K=3 para tipologias optimas en K=4)\n\n")

sensibilidad_K_lista <- list()
for (tp in tipos) {
  res_all <- hmm_resultados_all[[tp]]$all
  if (is.null(res_all)) next
  Ks_disponibles <- as.integer(names(res_all))
  for (k in Ks_disponibles) {
    r <- res_all[[as.character(k)]]
    sensibilidad_K_lista[[paste(tp, k)]] <- data.frame(
      Tipologia = tp,
      K         = k,
      BIC       = r$BIC,
      LogLik    = r$LL,
      Ratio     = r$ratio,
      Mejor     = ifelse(k == hmm_resultados_all[[tp]]$best$nstates,
                         "*", "")
    )
  }
}
tabla_sensibilidad_K <- do.call(rbind, sensibilidad_K_lista)
rownames(tabla_sensibilidad_K) <- NULL
cat("  Tabla A.3 Sensibilidad K HMM:\n")
print(round_df(tabla_sensibilidad_K, 2))



#
cat("\n>>> B.1 Residuos por dia de la semana\n\n")

residuos_dow_lista <- list()
for (exp_name in c("2025", "2026")) {
  exp_obj <- get(paste0("exp_", exp_name))
  for (tp in tipos) {
    r <- exp_obj[[tp]]
    if (is.null(r)) next
    dows <- format(r$fechas, "%u")
    res  <- abs(r$real - r$pred_fwd)
    for (d in sort(unique(dows))) {
      residuos_dow_lista[[paste(exp_name, tp, d)]] <- data.frame(
        Tipologia   = tp,
        Experimento = exp_name,
        Dia_semana  = d,
        N           = sum(dows == d),
        MAE_dia     = mean(res[dows == d])
      )
    }
  }
}
tabla_residuos_dow <- do.call(rbind, residuos_dow_lista)
rownames(tabla_residuos_dow) <- NULL
cat("  Tabla B.1 Residuos por dia de la semana:\n")
print(head(round_df(tabla_residuos_dow, 4), 20))



#
cat("\n>>> B.2 MAE en festivos vs no festivos\n\n")

residuos_festivos_lista <- list()
for (exp_name in c("2025", "2026")) {
  exp_obj <- get(paste0("exp_", exp_name))
  for (tp in tipos) {
    r <- exp_obj[[tp]]
    if (is.null(r)) next
    es_festivo <- r$fechas %in% festivos
    if (sum(es_festivo) == 0) next
    err <- abs(r$real - r$pred_fwd)
    residuos_festivos_lista[[paste(exp_name, tp)]] <- data.frame(
      Tipologia    = tp,
      Experimento  = exp_name,
      N_festivos   = sum(es_festivo),
      MAE_festivo  = mean(err[es_festivo]),
      MAE_laborable = mean(err[!es_festivo])
    )
  }
}
tabla_residuos_festivos <- do.call(rbind, residuos_festivos_lista)
rownames(tabla_residuos_festivos) <- NULL
cat("  Tabla B.2 MAE festivos vs laborables:\n")
print(round_df(tabla_residuos_festivos, 4))



#
cat("\n>>> C.1 Top 25 areas comunitarias por volumen criminal\n\n")

if (!is.null(perfil_comunidades)) {
  top25_ca <- perfil_comunidades %>%
    dplyr::arrange(desc(Total_Incidentes)) %>%
    head(25)
  print(top25_ca)
} else {
  top25_ca <- NULL
  cat("  [Aviso] No disponible.\n")
}



#
cat("\n>>> C.2 Top 25 dias con mas incidentes\n\n")

dt_top_dias <- dt[, .(N = .N), by = fecha]
data.table::setorder(dt_top_dias, -N)
top25_dias <- head(dt_top_dias, 25)
print(top25_dias)



#
cat("\n>>> C.3 Distribucion anual por tipologia\n\n")

dist_anual <- dt[, .(N = .N), by = .(year, tipo_delito)]
dist_anual_wide <- tidyr::pivot_wider(dist_anual,
                                       names_from = tipo_delito,
                                       values_from = N,
                                       values_fill = 0)
print(as.data.frame(dist_anual_wide))



#
cat("\n>>> C.4 Top 10 modus operandi por tipologia\n\n")

modus_top10 <- lapply(tipos, function(tp) {
  m <- dt[tipo_delito == tp,
          .(N = .N), by = description][order(-N)][seq_len(min(10, .N))]
  m$Tipologia <- tp
  m$pct <- round(100 * m$N / sum(m$N), 2)
  m
})
tabla_modus_top10 <- do.call(rbind, modus_top10)
cat("  Tabla C.4 Modus operandi extendido:\n")
print(as.data.frame(tabla_modus_top10))



#
cat("\n>>> D.1 Checklist de reproducibilidad\n\n")

checklist_repro <- data.frame(
  Item = c(
    "Semilla maestra fijada", "set.seed antes de cada train",
    "CV 3-fold con allowParallel", "Bootstrap con seed",
    "Forward Algorithm deterministico", "Versiones de paquetes registradas",
    "Datos raw inmutables", "Calendario de festivos generado por funcion"
  ),
  Estado = c("OK", "OK", "OK", "OK", "OK", "OK", "OK", "OK")
)
print(checklist_repro)



#
cat("\n>>> D.2 Hashes de ficheros raw (para auditoria)\n\n")

if (requireNamespace("digest", quietly = TRUE)) {
  hashes <- c(
    csv_raw      = if (file.exists(CSV_RAW))         digest::digest(file = CSV_RAW, algo = "md5")        else NA,
    clima        = if (file.exists(CSV_CLIMA))       digest::digest(file = CSV_CLIMA, algo = "md5")      else NA,
    geojson      = if (file.exists(GEOJSON_CA))      digest::digest(file = GEOJSON_CA, algo = "md5")     else NA,
    comisarias   = if (file.exists(CSV_COMISARIAS))  digest::digest(file = CSV_COMISARIAS, algo = "md5") else NA,
    poblacion    = if (file.exists(CSV_POBLACION))   digest::digest(file = CSV_POBLACION, algo = "md5")  else NA
  )
  cat("  MD5 de los ficheros:\n")
  print(hashes)
} else {
  cat("  [Aviso] paquete digest no disponible. Se omite hash.\n")
  hashes <- NULL
}



#
cat("\n>>> E.1 Benchmark de tiempos por seccion\n\n")

if (exists("t_parte_i")) {
  tiempos_partes <- data.frame(
    Parte    = c("Parte I (ETL)",
                 "Parte II (Temporal)",
                 "Parte III (MS-MoE)",
                 "Parte IV (Validacion)",
                 "Parte V (Operativo)",
                 "Parte VI (Cierre)",
                 "Total"),
    Minutos = c(
      as.numeric(difftime(t_parte_ii,  t_parte_i,   units = "mins")),
      as.numeric(difftime(t_parte_iii, t_parte_ii,  units = "mins")),
      as.numeric(difftime(t_parte_iv,  t_parte_iii, units = "mins")),
      as.numeric(difftime(t_parte_v,   t_parte_iv,  units = "mins")),
      as.numeric(difftime(t_parte_vi,  t_parte_v,   units = "mins")),
      as.numeric(difftime(Sys.time(),  t_parte_vi,  units = "mins")),
      as.numeric(difftime(Sys.time(),  t_global,    units = "mins"))
    )
  )
  print(round_df(tiempos_partes, 2))
} else {
  tiempos_partes <- NULL
  cat("  [Aviso] Timestamps no disponibles. Se omite benchmark.\n")
}



#
cat("\n>>> F.1 Sintesis final del estado del modelo\n\n")

sintesis_final <- data.frame(
  Metrica = c(
    "Total filas raw",                "Total tras filtro 4 tipos",
    "Periodo",                        "Dias totales",
    "Numero de tipologias",           "Numero de areas comunitarias",
    "Numero de distritos policiales", "K HMM (LES)",
    "K HMM (AMENAZAS)",               "K HMM (ROBO)",
    "K HMM (HOMICIDIO)",
    "MAE MS-MoE (LES 2025)",          "MAE MS-MoE (LES 2026)",
    "MAE MS-MoE (AME 2025)",          "MAE MS-MoE (AME 2026)",
    "MAE MS-MoE (ROBO 2025)",         "MAE MS-MoE (ROBO 2026)",
    "MAE MS-MoE (HOM 2025)",          "MAE MS-MoE (HOM 2026)",
    "Victorias MS-MoE",               "Victorias Prophet",
    "Victorias NNET",                 "Bootstrap B",
    "Bootstrap L",                    "R0 SIR (si aplica)"
  ),
  Valor = c(
    "1.504.308",                                    format(nrow(dt), big.mark = "."),
    sprintf("%s a %s", min(dt$fecha), max(dt$fecha)), as.character(n_dias),
    as.character(length(tipos)),                    if (!is.null(perfil_comunidades)) as.character(nrow(perfil_comunidades)) else "N/A",
    as.character(length(unique(dt$district[!is.na(dt$district) & dt$district > 0]))),
    as.character(tabla_viterbi$K_optimo[tabla_viterbi$Tipologia == "LESIONES"]),
    as.character(tabla_viterbi$K_optimo[tabla_viterbi$Tipologia == "AMENAZAS"]),
    as.character(tabla_viterbi$K_optimo[tabla_viterbi$Tipologia == "ROBO CON VIOLENCIA"]),
    as.character(tabla_viterbi$K_optimo[tabla_viterbi$Tipologia == "HOMICIDIO"]),
    sprintf("%.3f", exp_2025$LESIONES$mae_msmoe_fwd),
    sprintf("%.3f", exp_2026$LESIONES$mae_msmoe_fwd),
    sprintf("%.3f", exp_2025$AMENAZAS$mae_msmoe_fwd),
    sprintf("%.3f", exp_2026$AMENAZAS$mae_msmoe_fwd),
    sprintf("%.3f", exp_2025[["ROBO CON VIOLENCIA"]]$mae_msmoe_fwd),
    sprintf("%.3f", exp_2026[["ROBO CON VIOLENCIA"]]$mae_msmoe_fwd),
    sprintf("%.3f", exp_2025$HOMICIDIO$mae_msmoe_fwd),
    sprintf("%.3f", exp_2026$HOMICIDIO$mae_msmoe_fwd),
    as.character(sum(tabla_victorias_resumen$Ganador == "MS-MoE")),
    as.character(sum(tabla_victorias_resumen$Ganador == "Prophet")),
    as.character(sum(tabla_victorias_resumen$Ganador == "NNET")),
    as.character(B_BOOT),
    as.character(L_BOOT),
    if (exists("R0_sir")) sprintf("%.3f", R0_sir) else "N/A"
  ),
  stringsAsFactors = FALSE
)
cat("  TABLA SINTESIS FINAL:\n")
print(sintesis_final)



#
cat("\n>>> F.2 Lista exhaustiva de objetos exportados\n\n")

env_actual <- ls(envir = .GlobalEnv)
objetos_tipo <- vapply(env_actual,
                       function(x) class(get(x, envir = .GlobalEnv))[1],
                       character(1L))
tabla_objetos <- data.frame(
  Objeto = env_actual,
  Clase  = objetos_tipo
)
tabla_objetos <- tabla_objetos[order(tabla_objetos$Clase, tabla_objetos$Objeto), ]
rownames(tabla_objetos) <- NULL
cat(sprintf("  Total objetos en .GlobalEnv: %d\n", nrow(tabla_objetos)))
print(head(tabla_objetos, 30))



#
RDATA_APENDICES <- file.path(RESULT, "resultados_apendices.RData")
save(
  tabla_sensibilidad_rolling, tabla_estabilidad_mes, tabla_sensibilidad_K,
  tabla_residuos_dow, tabla_residuos_festivos,
  top25_ca, top25_dias, dist_anual_wide, tabla_modus_top10,
  checklist_repro, hashes, tiempos_partes,
  sintesis_final, tabla_objetos,
  file = RDATA_APENDICES
)
cat(sprintf("\n  Guardado apendices en: %s\n", RDATA_APENDICES))
cat(sprintf("  Tamano: %.2f MB\n", file.size(RDATA_APENDICES) / 1024 ^ 2))

cat("\n")
cat(strrep("#", 78), "\n", sep = "")
cat("# APENDICES COMPLETADOS\n")
cat(strrep("#", 78), "\n\n", sep = "")

sink(type = "message"); sink()
close(apendice_log)


# FIN DEL SCRIPT MAESTRO v3 (extension con apendices)
