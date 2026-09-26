# HistoWUI raster mapping scripts

These scripts create historical and contemporary HistoWUI rasters from prepared building-count and NLCD inputs. `main.py` connects the processing steps and keeps the historical 2020 output separate from the contemporary 2020 output. The repository contains the mapping workflow; it does **not** create the input building inventories, assign construction years, prepare NLCD, or split and align the input rasters.

## Workflow

| Stage | Script | Result |
| --- | --- | --- |
| Neighborhood calculation | `Step1_Function_density_coverage.py`, called by `Step3_building_den.py` and `Step3_vegetation_cov.py` | Local building counts and wildland vegetation coverage |
| Vegetation patch and buffer | `Step2_identify_wv_cluster_and 24buffer_for_interfaceWUI.py` | Vegetation patches of at least 5 km² and their 2.4 km buffer |
| WUI classification | `Step3_identify_WUI.py` | Intermix, interface, and other land classes |
| Run the complete workflow | `main.py` | Configures paths, filenames, years, and processing order |

The script names reflect their development history. Use `main.py` to run them in the correct order.

## Requirements

- Python 3.9 or newer with NumPy, SciPy, and Rasterio (see `requirements.txt`).
- GDAL command-line programs `gdalbuildvrt` and `gdal_translate` on `PATH`.
- Enough disk space for the intermediate GeoTIFFs and enough RAM to hold a full annual CONUS mosaic during patch filtering and WUI classification.

One possible installation is:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
gdalbuildvrt --version
gdal_translate --version
```

Install the GDAL command-line programs through your system package manager or your computing environment if they are unavailable. Rasterio installation can also depend on the GDAL libraries available on that system.

## Prepare the inputs

All rasters must be single-band, 30 m rasters. The building and NLCD parts must use compatible pixel grids and have a spatial overlap of at least 19 pixels at internal part boundaries. For the final WUI step, the mosaicked building, vegetation-coverage, vegetation-buffer, and full-extent NLCD rasters must have the same CRS, shape, transform, and extent. The script checks this alignment before classification; it does not align the inputs for you.

`main.py` takes three input directories:

1. **Building parts:** per-pixel building counts. Historical filenames must contain a four-digit year and end in `_part1.tif`, `_part2.tif`, or `_part3.tif`. Contemporary filenames must additionally start with `contemporary_`. For example, `BuildingCountPerPixel_ALL_1990_30m_part1.tif` and `contemporary_BuildingCountPerPixel_ALL_2020_30m_part1.tif`.
2. **NLCD parts:** annual NLCD class rasters whose filenames contain a four-digit year and exactly one `part1`, `part2`, or `part3` identifier. For example, `CONUS_NLCD_1990_part1.tif`.
3. **NLCD reference rasters:** one full-extent raster per year, named `CONUS_NLCD_{year}_reextent.tif` by default. For example, `CONUS_NLCD_1990_reextent.tif`. Use `--nlcd-reference-pattern` if your names differ.

The code expects the three part identifiers `part1`, `part2`, and `part3` for every selected year. The historical building inventory is the subset used for the time series; the `contemporary_` inventory is used for the separate complete 2020 map. Supply the input data and document their sources and permissions separately when releasing the repository.

## Run

From the repository directory, preview the planned steps without reading or writing data:

```bash
python main.py \
  --building-parts /data/building_count_per_pixel/ALL_split \
  --nlcd-parts /data/NLCD_conus_split \
  --nlcd-reference /data/NLCD_conus_reextent \
  --output-root /data/HistoWUI_output \
  --mode both \
  --dry-run
```

Remove `--dry-run` to execute the workflow. `--mode historical` is the default and processes 1985, 1990, …, 2020. `--mode contemporary` creates only the complete 2020 map. `--mode both` creates both series and reuses the same 2020 vegetation inputs. To process selected historical years, pass `--years 1990,2000,2010,2020`. The `--years` option affects historical processing only; the contemporary inventory always uses 2020. Existing intermediate and final outputs are skipped unless `--overwrite` is supplied.

For example, to make only the historical 1990 map:

```bash
python main.py \
  --building-parts /data/building_count_per_pixel/ALL_split \
  --nlcd-parts /data/NLCD_conus_split \
  --nlcd-reference /data/NLCD_conus_reextent \
  --output-root /data/HistoWUI_output \
  --mode historical \
  --years 1990
```

## Mapping settings

The published entry point uses a circular radius of **19 pixels**, or **570 m** on the 30 m grid. The discrete kernel has 1,129 pixels, equivalent to 1.0161 km². The building output is a neighborhood **count**, although its legacy filenames contain `bdgden`. The NLCD neighborhood output is scaled from 0 to 10,000.

The WUI classification retains the thresholds in `Step3_identify_WUI.py`: building count greater than 6.17, high building count at least 49.42, and vegetation coverage greater than 5,000 for intermix WUI. Interface WUI uses the same building threshold, vegetation coverage at most 5,000, and the vegetation patch buffer. The buffer step retains 8-connected patches with at least 5 km² of vegetation coverage at or above 75% and includes cells within 2.4 km of those patches. The final categorical raster is reprojected to EPSG:5070 at 30 m using nearest-neighbor resampling.

| Value | Class |
| ---: | --- |
| 0 | NoData or unclassified |
| 1 | Intermix WUI |
| 2 | Interface WUI |
| 3 | Vegetated, low building count |
| 4 | Vegetated, no buildings |
| 5 | Non-vegetated, low building count |
| 6 | Non-vegetated, high building count |
| 7 | Water or perennial snow/ice |

Keep the radius and classification thresholds together when designing a sensitivity analysis. `main.py` fixes them to the baseline settings so a radius change cannot silently reuse the baseline building thresholds.

## Outputs

The output root contains:

```text
building_density_570m_conus_split/ALL/        local building-count parts
Vegetation_coverage_570m_conus_split/         local vegetation-coverage parts
vegetation_buffer/                            annual binary mosaics and buffers
HistoWUI_time_series/HistoWUI_ALL_<year>.tif  historical WUI rasters
HistoWUI_contemporary_2020/
    HistoWUI_ALL_2020_compl.tif                complete 2020 WUI raster
```

The complete 2020 output has a distinct filename and directory because the historical and contemporary 2020 building inventories differ.

## Data handling and reproducibility

The neighborhood function treats building NoData as zero when summing counts; unknown NLCD class values do not contribute to wildland vegetation coverage. Check the valid-data footprints of prepared inputs before comparing areas or fractions. The final class value `0` also includes cells left unclassified by the criteria, so it should not automatically be interpreted as a measured non-WUI class.

The patch and WUI scripts load full annual rasters into memory. For CONUS runs, use a compute node with sufficient memory and local scratch space. Record the input dataset versions, preprocessing steps, and software versions used for the manuscript analysis alongside the repository release.
