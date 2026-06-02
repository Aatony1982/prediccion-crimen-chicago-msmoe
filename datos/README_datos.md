# Datos

Este directorio contiene los ficheros auxiliares de menor tamaño, ya incluidos
en el repositorio. El **histórico delictivo completo no se versiona** por su
tamaño (supera el límite de GitHub) y debe descargarse aparte.

## Ficheros incluidos

| Fichero | Descripción |
|---------|-------------|
| `clima_diario.csv` | Temperatura media diaria (Open-Meteo, estación O'Hare). |
| `comisarias_cpd.csv` | Coordenadas y plantilla de las 22 comisarías del CPD. |
| `community_areas.geojson` | Polígonos de las 77 áreas comunitarias de Chicago. |
| `chicago_population.csv` | Población, renta y composición por área comunitaria. |

## Fichero que debes descargar

`delitos_raw_download_socrata.csv` — histórico de incidentes del Chicago Police
Department.

### Pasos

1. Abre el dataset **Crimes - 2001 to Present** del Chicago Data Portal:
   <https://data.cityofchicago.org/Public-Safety/Crimes-2001-to-Present/ijzp-q8t2>
   (identificador Socrata `ijzp-q8t2`).

2. Exporta a **CSV**. Puedes descargar el histórico completo o, para reducir el
   tamaño, filtrar por fecha al periodo **2020-01-01 a 2026-05-31** antes de
   exportar. El script trabaja internamente con ese rango.

   Alternativa por API (sin navegador):

   ```bash
   curl "https://data.cityofchicago.org/resource/ijzp-q8t2.csv?$where=date between '2020-01-01T00:00:00' and '2026-05-31T23:59:59'&$limit=2000000" -o delitos_raw_download_socrata.csv
   ```

3. Guarda el fichero en esta carpeta con el nombre exacto:

   ```
   datos/delitos_raw_download_socrata.csv
   ```

El script espera, como mínimo, las columnas `case_number`, `date`,
`primary_type`, `year`, `latitude`, `longitude`, `community_area`, `district`,
`arrest`, `domestic`, `description` y `location_description`.
