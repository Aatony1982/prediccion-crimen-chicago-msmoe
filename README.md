# Predicción de crimen urbano en Chicago — MS-MoE

Código y datos para **reproducir íntegramente** el experimento de predicción
diaria de criminalidad violenta en Chicago mediante un modelo
**Markov-Switching Mixture of Experts (MS-MoE)**: un Modelo Oculto de Markov
(HMM) gaussiano con selección de estados por BIC, dos expertos de red neuronal
promediada (avNNet) especializados por régimen, y un selector basado en el
*Forward Algorithm* de Hamilton.

El modelo se entrena y valida sobre cuatro tipologías violentas (lesiones,
amenazas, robo con violencia y homicidio) en dos periodos de test
independientes (2025 y 2026), comparándolo frente a Prophet y avNNet, y se
extiende a una capa operativa (grafo de comisarías, Dijkstra, asignación de
patrullas por turno).

## Estructura del repositorio

```
.
├── reproducir_experimento.R     Script único que reproduce todo el análisis
├── datos/
│   ├── clima_diario.csv          Temperatura diaria (Open-Meteo)
│   ├── comisarias_cpd.csv        Coordenadas de las 22 comisarías del CPD
│   ├── community_areas.geojson   Límites de las 77 áreas comunitarias
│   ├── chicago_population.csv    Población y renta por área comunitaria
│   └── README_datos.md           Cómo descargar el histórico delictivo completo
├── resultados/                   Salida generada por el script (no versionada)
├── AVISO_DERECHOS.txt            Aviso de derechos de autor (todos reservados)
└── README.md
```

## Requisitos

- **R 4.4** o superior.
- Los paquetes necesarios se **instalan automáticamente** la primera vez
  (data.table, dplyr, tidyr, lubridate, caret, nnet, forecast, tseries,
  changepoint, depmixS4, MASS, pscl, doParallel, foreach, sf, spdep, geosphere,
  igraph, ggplot2, zoo, prophet, digest, osmdata, entre otros).
- Hardware recomendado: 7+ núcleos lógicos y 16 GB de RAM.

## Cómo reproducir

1. **Descargar el histórico delictivo** y colocarlo en `datos/`
   (ver instrucciones en [`datos/README_datos.md`](datos/README_datos.md)).
   El resto de ficheros de datos ya están incluidos.

2. **Ejecutar el script** desde la raíz del repositorio:

   ```bash
   Rscript reproducir_experimento.R
   ```

   El script usa rutas relativas y emplea una **semilla aleatoria fija (2026)**,
   de modo que ejecuciones sucesivas reproducen resultados idénticos. Tiempo
   aproximado: 30–50 minutos.

3. **Salida**: `resultados/resultados.RData`, que contiene todos los objetos
   numéricos, tablas y figuras del análisis.

## Fuentes de datos

| Dato | Fuente |
|------|--------|
| Histórico delictivo | [Chicago Data Portal — Crimes (Socrata, dataset `ijzp-q8t2`)](https://data.cityofchicago.org/Public-Safety/Crimes-2001-to-Present/ijzp-q8t2) |
| Temperatura diaria | [Open-Meteo](https://open-meteo.com/) |
| Áreas comunitarias | [Chicago Data Portal — Community Areas](https://data.cityofchicago.org/) |
| Población y renta por área comunitaria | Censo (American Community Survey) |

> La población y la renta **no** se utilizan como predictores del modelo
> (diseño deliberadamente ciego a la demografía). Se incluyen únicamente para
> la verificación de equidad del análisis.

## Derechos de autor

**Todos los derechos reservados** (ver [`AVISO_DERECHOS.txt`](AVISO_DERECHOS.txt)).
Este repositorio se publica únicamente con fines de transparencia y
reproducibilidad académica. No se concede licencia de uso, copia, modificación ni
redistribución del código sin autorización expresa y por escrito del autor. Los
datos pertenecen a sus respectivas fuentes y se rigen por sus propias condiciones
de uso.
