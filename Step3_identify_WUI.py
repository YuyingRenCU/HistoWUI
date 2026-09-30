#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create annual HistoWUI classification rasters from prepared input layers.

The script combines building-density, vegetation-coverage, vegetation-buffer,
and NLCD rasters for each requested year. Input rasters must already use the
same grid, extent, resolution, and coordinate reference system. The classified
output is reprojected to EPSG:5070 at 30 m using nearest-neighbor resampling.

Class values
------------
0 : NoData or unclassified
1 : Intermix WUI
2 : Interface WUI
3 : Wildland - Low Building
4 : Wildland - No Building
5 : Other Vegetated (Urban, Agriculture, Barren)
6 : Urban - High Building
7 : Water or perennial snow/ice

Important
---------
The building-density thresholds below apply to the neighborhood used to
generate the input building-density rasters. Recalculate the thresholds if the
neighborhood definition changes.

Class names follow the published data legend. Classes 5 and 6 are assigned
from building counts and local wildland-vegetation coverage; they do not
require a separate NLCD urban, agriculture, or barren class filter.
"""

from pathlib import Path

import numpy as np
import rasterio as rio
from rasterio.crs import CRS
from rasterio.merge import merge
from rasterio.transform import array_bounds
from rasterio.warp import Resampling, calculate_default_transform, reproject


# -----------------------------------------------------------------------------
# User configuration
# -----------------------------------------------------------------------------
# Replace these example paths with directories containing the actual inputs.
BUILDING_TYPE = "ALL"

VEGETATION_COVERAGE_DIR = Path("/path/to/vegetation_coverage_rasters")
VEGETATION_BUFFER_DIR = Path("/path/to/vegetation_buffer_rasters")
BUILDING_DENSITY_DIR = Path("/path/to/building_density_rasters") / BUILDING_TYPE
NLCD_DIR = Path("/path/to/nlcd_rasters")
OUTPUT_DIR = Path("/path/to/histowui_output")

# Update these patterns if the input filenames use a different convention.
# Each building-density and vegetation-coverage year is expected to contain
# three parts named part1, part2, and part3.
BUILDING_DENSITY_PATTERN = "conus_*{year}*{part}.tif"
VEGETATION_COVERAGE_PATTERN = "*{year}*{part}.tif"
VEGETATION_BUFFER_PATTERN = "*{year}*buffer*{distance}*.tif"
NLCD_PATTERN = "CONUS_NLCD_{year}_reextent.tif"

YEARS = [str(year) for year in range(1985, 2025, 5)]
PARTS = ("part1", "part2", "part3")


# -----------------------------------------------------------------------------
# Classification parameters
# -----------------------------------------------------------------------------
BUILDING_DENSITY_THRESHOLD = 6.17
HIGH_BUILDING_DENSITY_THRESHOLD = 49.42
VEGETATION_COVERAGE_THRESHOLD = 5000

CONTINUOUS_RASTER_NODATA = -999
NLCD_NODATA = 250
VEGETATION_BUFFER_VALUE = 1
VEGETATION_BUFFER_DISTANCE_M = 2400

WATER_AND_ICE_CLASSES = (11, 12)

TARGET_CRS = CRS.from_epsg(5070)
TARGET_RESOLUTION_M = 30
OUTPUT_NODATA = 0


def find_one(directory: Path, pattern: str) -> Path:
    """Return the only file matching a pattern and report ambiguous inputs."""
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one file matching '{pattern}' in '{directory}', "
            f"but found {len(matches)}: {matches}"
        )
    return matches[0]


def merge_year_parts(directory: Path, pattern: str, year: str):
    """Merge the three raster parts for one year using the maximum value."""
    paths = [
        find_one(directory, pattern.format(year=year, part=part))
        for part in PARTS
    ]
    sources = [rio.open(path) for path in paths]
    try:
        mosaic, transform = merge(sources, method="max")
        crs = sources[0].crs
    finally:
        for source in sources:
            source.close()

    return mosaic[0], transform, crs


def read_single_band(path: Path):
    """Read a single-band raster and return its array, profile, and transform."""
    with rio.open(path) as source:
        array = source.read(1)
        profile = source.profile.copy()
        transform = source.transform
        crs = source.crs
    return array, profile, transform, crs


def validate_common_grid(reference_array, reference_transform, reference_crs, layers):
    """Verify that all input arrays use the reference raster grid."""
    for name, array, transform, crs in layers:
        if array.shape != reference_array.shape:
            raise ValueError(
                f"{name} shape {array.shape} does not match the NLCD shape "
                f"{reference_array.shape}"
            )
        if not transform.almost_equals(reference_transform):
            raise ValueError(f"{name} transform does not match the NLCD transform")
        if crs != reference_crs:
            raise ValueError(f"{name} CRS {crs} does not match the NLCD CRS {reference_crs}")


def classify_wui(building_density, vegetation_coverage, vegetation_buffer, nlcd):
    """Apply the HistoWUI classification criteria without changing their order."""
    classified = np.zeros(nlcd.shape, dtype=np.uint8)

    vegetated = vegetation_coverage > VEGETATION_COVERAGE_THRESHOLD
    nonvegetated_valid = (
        (vegetation_coverage <= VEGETATION_COVERAGE_THRESHOLD)
        & (vegetation_coverage > CONTINUOUS_RASTER_NODATA)
    )

    # Non-WUI classes
    classified[(building_density == 0) & vegetated] = 4
    classified[
        (building_density > 0)
        & (building_density <= BUILDING_DENSITY_THRESHOLD)
        & vegetated
    ] = 3
    classified[
        (building_density >= 0)
        & (building_density < HIGH_BUILDING_DENSITY_THRESHOLD)
        & nonvegetated_valid
        & (nlcd != NLCD_NODATA)
    ] = 5
    classified[
        (building_density >= HIGH_BUILDING_DENSITY_THRESHOLD)
        & nonvegetated_valid
    ] = 6

    # WUI classes
    classified[
        vegetated & (building_density > BUILDING_DENSITY_THRESHOLD)
    ] = 1
    classified[
        (building_density > BUILDING_DENSITY_THRESHOLD)
        & nonvegetated_valid
        & (vegetation_buffer == VEGETATION_BUFFER_VALUE)
    ] = 2

    # Water and perennial snow/ice override all preceding classes.
    classified[np.isin(nlcd, WATER_AND_ICE_CLASSES)] = 7

    # Set up the nodata extent for the final output rasters
    classified[(nlcd != NLCD_NODATA)] = 0

    return classified


def write_reprojected_output(classified, source_profile, output_path: Path):
    """Reproject a classified raster to EPSG:5070 and save it as a GeoTIFF."""
    source_crs = source_profile["crs"]
    source_transform = source_profile["transform"]
    source_height, source_width = classified.shape

    left, bottom, right, top = array_bounds(
        source_height,
        source_width,
        source_transform,
    )
    destination_transform, destination_width, destination_height = (
        calculate_default_transform(
            source_crs,
            TARGET_CRS,
            source_width,
            source_height,
            left,
            bottom,
            right,
            top,
            resolution=TARGET_RESOLUTION_M,
        )
    )

    output_profile = source_profile.copy()
    output_profile.update(
        driver="GTiff",
        crs=TARGET_CRS,
        transform=destination_transform,
        width=destination_width,
        height=destination_height,
        count=1,
        dtype=rio.uint8,
        nodata=OUTPUT_NODATA,
        compress="LZW",
        tiled=True,
        BIGTIFF="YES",
    )

    with rio.open(output_path, "w", **output_profile) as destination:
        reproject(
            source=classified,
            destination=rio.band(destination, 1),
            src_transform=source_transform,
            src_crs=source_crs,
            dst_transform=destination_transform,
            dst_crs=TARGET_CRS,
            src_nodata=OUTPUT_NODATA,
            dst_nodata=OUTPUT_NODATA,
            resampling=Resampling.nearest,
        )


def process_year(year: str):
    """Read, validate, classify, and save the inputs for one year."""
    print(f"Processing {year}", flush=True)

    building_density, building_transform, building_crs = merge_year_parts(
        BUILDING_DENSITY_DIR,
        BUILDING_DENSITY_PATTERN,
        year,
    )
    vegetation_coverage, vegetation_transform, vegetation_crs = merge_year_parts(
        VEGETATION_COVERAGE_DIR,
        VEGETATION_COVERAGE_PATTERN,
        year,
    )

    buffer_path = find_one(
        VEGETATION_BUFFER_DIR,
        VEGETATION_BUFFER_PATTERN.format(
            year=year,
            distance=VEGETATION_BUFFER_DISTANCE_M,
        ),
    )
    vegetation_buffer, _, buffer_transform, buffer_crs = read_single_band(
        buffer_path
    )

    nlcd_path = find_one(NLCD_DIR, NLCD_PATTERN.format(year=year))
    nlcd, nlcd_profile, nlcd_transform, nlcd_crs = read_single_band(nlcd_path)

    validate_common_grid(
        nlcd,
        nlcd_transform,
        nlcd_crs,
        [
            (
                "Building-density raster",
                building_density,
                building_transform,
                building_crs,
            ),
            (
                "Vegetation-coverage raster",
                vegetation_coverage,
                vegetation_transform,
                vegetation_crs,
            ),
            (
                "Vegetation-buffer raster",
                vegetation_buffer,
                buffer_transform,
                buffer_crs,
            ),
        ],
    )

    classified = classify_wui(
        building_density,
        vegetation_coverage,
        vegetation_buffer,
        nlcd,
    )
    values, counts = np.unique(classified, return_counts=True)
    print(
        "Class counts: "
        + ", ".join(f"{value}={count:,}" for value, count in zip(values, counts)),
        flush=True,
    )

    output_path = OUTPUT_DIR / f"HistoWUI_{BUILDING_TYPE}_{year}.tif"
    write_reprojected_output(classified, nlcd_profile, output_path)
    print(f"Saved {output_path}", flush=True)


def main():
    """Create HistoWUI rasters for all configured years."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for year in YEARS:
        process_year(year)


if __name__ == "__main__":
    main()
