#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create buffered high-vegetation patches from annual coverage rasters.

For each year, the workflow:

1. Reclassifies vegetation coverage to a binary raster using a 75% threshold.
2. Mosaics the input chunks.
3. Identifies eight-neighbor-connected vegetation patches of at least 5 km².
4. Creates a 2.4 km Euclidean buffer around the retained patches.
5. Writes either the vegetation core plus buffer or the external buffer ring.

Input rasters must be single-band, aligned 30 m rasters. Vegetation coverage is
assumed to be stored on a 0--10,000 scale, so 7,500 represents 75%. The script
requires GDAL command-line programs ``gdalbuildvrt`` and ``gdal_translate``.
"""

import gc
import shutil
import subprocess
from pathlib import Path
from time import perf_counter

import numpy as np
import rasterio as rio
from scipy.ndimage import distance_transform_edt, generate_binary_structure, label


# -----------------------------------------------------------------------------
# User configuration
# -----------------------------------------------------------------------------
# Replace these example paths with the actual input and output directories.
INPUT_CHUNK_DIR = Path("/path/to/annual_vegetation_coverage_chunks")
OUTPUT_ROOT = Path("/path/to/vegetation_buffer_output")

# Change this pattern if the input filenames use a different convention.
# It must match every chunk for the requested year.
INPUT_CHUNK_PATTERN = "vegetation_coverage_{year}_part*.tif"
EXPECTED_CHUNKS = 3

# Process one year or every year in the inclusive range below.
RUN_ONLY_ONE_YEAR = True
TARGET_YEAR = 2020
START_YEAR = 1986
END_YEAR = 2023

# True writes retained vegetation patches and their buffer. False writes only
# the external buffer ring.
INCLUDE_VEGETATED_CORE = True

# Existing outputs are retained unless their corresponding option is True.
OVERWRITE_BINARY_CHUNKS = False
OVERWRITE_MOSAIC = False
OVERWRITE_FINAL_OUTPUT = False


# -----------------------------------------------------------------------------
# Output directories
# -----------------------------------------------------------------------------
BINARY_CHUNK_DIR = OUTPUT_ROOT / "binary_wv75_chunks"
BINARY_MOSAIC_DIR = OUTPUT_ROOT / "binary_wv75_annual_mosaics"
FINAL_OUTPUT_DIR = OUTPUT_ROOT / "wv75_patch5km2_buffer2400m"


# -----------------------------------------------------------------------------
# Analysis parameters
# -----------------------------------------------------------------------------
VEGETATION_COVERAGE_THRESHOLD = 7500
PIXEL_SIZE_M = 30
MIN_PATCH_AREA_M2 = 5_000_000
MIN_PATCH_PIXELS = int(np.ceil(MIN_PATCH_AREA_M2 / PIXEL_SIZE_M**2))
BUFFER_DISTANCE_M = 2400

OUTPUT_NODATA = 255
OUTPUT_BLOCK_SIZE = 512


def nodata_is_nan(nodata):
    """Return True when a raster NoData value is NaN."""
    return (
        nodata is not None
        and isinstance(nodata, (float, np.floating))
        and np.isnan(nodata)
    )


def run_command(command):
    """Run an external command and stop if it fails."""
    command = [str(value) for value in command]
    print("Running:", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def require_gdal_commands():
    """Confirm that the required GDAL programs are available."""
    missing = [
        command
        for command in ("gdalbuildvrt", "gdal_translate")
        if shutil.which(command) is None
    ]
    if missing:
        raise RuntimeError(
            "Required GDAL command(s) not found on PATH: " + ", ".join(missing)
        )


def years_to_process():
    """Return the configured processing years as strings."""
    if RUN_ONLY_ONE_YEAR:
        return [str(TARGET_YEAR)]
    return [str(year) for year in range(START_YEAR, END_YEAR + 1)]


def grid_origins_align(reference_transform, transform, resolution):
    """Check whether two raster origins fall on the same pixel grid."""
    x_offset = (transform.c - reference_transform.c) / resolution[0]
    y_offset = (transform.f - reference_transform.f) / abs(resolution[1])
    return np.isclose(x_offset, round(x_offset)) and np.isclose(
        y_offset, round(y_offset)
    )


def inspect_and_validate_chunks(chunk_paths):
    """Validate band count, CRS, resolution, and alignment of input chunks."""
    metadata = []
    for path in chunk_paths:
        with rio.open(path) as source:
            metadata.append(
                {
                    "path": path,
                    "crs": source.crs,
                    "res": source.res,
                    "dtype": source.dtypes[0],
                    "count": source.count,
                    "nodata": source.nodata,
                    "transform": source.transform,
                    "bounds": source.bounds,
                }
            )

    reference = metadata[0]
    expected_resolution = (PIXEL_SIZE_M, PIXEL_SIZE_M)

    for item in metadata:
        if item["count"] != 1:
            raise ValueError(
                f"Expected one raster band in {item['path']}, "
                f"but found {item['count']}"
            )
        if item["crs"] is None:
            raise ValueError(f"Input raster has no CRS: {item['path']}")
        if not np.allclose(item["res"], expected_resolution):
            raise ValueError(
                f"Expected {expected_resolution} m pixels in {item['path']}, "
                f"but found {item['res']}"
            )
        if item["crs"] != reference["crs"]:
            raise ValueError(
                f"CRS mismatch between {reference['path']} and {item['path']}"
            )
        if not np.allclose(item["res"], reference["res"]):
            raise ValueError(
                f"Resolution mismatch between {reference['path']} and "
                f"{item['path']}"
            )
        if not grid_origins_align(
            reference["transform"], item["transform"], reference["res"]
        ):
            raise ValueError(
                f"Pixel-grid mismatch between {reference['path']} and "
                f"{item['path']}"
            )

    print("Input chunks:")
    for item in metadata:
        print(
            f"  {item['path']}\n"
            f"    CRS: {item['crs']}\n"
            f"    resolution: {item['res']}\n"
            f"    dtype: {item['dtype']}\n"
            f"    NoData: {item['nodata']}\n"
            f"    bounds: {item['bounds']}"
        )


def reclassify_chunk_to_uint8(input_path, output_path, coverage_threshold):
    """Reclassify one coverage chunk to 1, 0, and NoData block by block."""
    if output_path.exists() and not OVERWRITE_BINARY_CHUNKS:
        print(f"Binary chunk exists; skipping: {output_path}")
        return

    print(f"Reclassifying {input_path}")
    start = perf_counter()

    with rio.open(input_path) as source:
        profile = source.profile.copy()
        input_nodata = source.nodata
        profile.update(
            dtype=rio.uint8,
            count=1,
            nodata=OUTPUT_NODATA,
            compress="deflate",
            zlevel=6,
            predictor=1,
            tiled=True,
            blockxsize=OUTPUT_BLOCK_SIZE,
            blockysize=OUTPUT_BLOCK_SIZE,
            BIGTIFF="IF_SAFER",
        )

        valid_total = 0
        high_coverage_total = 0

        with rio.open(output_path, "w", **profile) as destination:
            for _, window in source.block_windows(1):
                coverage = source.read(1, window=window)

                if input_nodata is None:
                    if np.issubdtype(coverage.dtype, np.floating):
                        valid = ~np.isnan(coverage)
                    else:
                        valid = np.ones(coverage.shape, dtype=bool)
                elif nodata_is_nan(input_nodata):
                    valid = ~np.isnan(coverage)
                else:
                    valid = coverage != input_nodata
                    if np.issubdtype(coverage.dtype, np.floating):
                        valid &= ~np.isnan(coverage)

                high_coverage = valid & (coverage >= coverage_threshold)
                binary = np.full(coverage.shape, OUTPUT_NODATA, dtype=np.uint8)
                binary[valid] = 0
                binary[high_coverage] = 1

                destination.write(binary, 1, window=window)
                valid_total += int(valid.sum())
                high_coverage_total += int(high_coverage.sum())

    percentage = (
        high_coverage_total / valid_total * 100 if valid_total > 0 else 0.0
    )
    print(
        f"  valid pixels: {valid_total:,}\n"
        f"  pixels meeting threshold: {high_coverage_total:,}\n"
        f"  percent meeting threshold: {percentage:.2f}%\n"
        f"  elapsed time: {(perf_counter() - start) / 60:.1f} minutes"
    )


def build_binary_mosaic(binary_chunk_paths, vrt_path, mosaic_path):
    """Create an annual VRT and compressed binary GeoTIFF mosaic."""
    if mosaic_path.exists() and not OVERWRITE_MOSAIC:
        print(f"Annual mosaic exists; skipping: {mosaic_path}")
        return

    vrt_path.unlink(missing_ok=True)
    run_command(
        [
            "gdalbuildvrt",
            "-overwrite",
            "-srcnodata",
            OUTPUT_NODATA,
            "-vrtnodata",
            OUTPUT_NODATA,
            "-resolution",
            "highest",
            vrt_path,
            *binary_chunk_paths,
        ]
    )

    if mosaic_path.exists() and OVERWRITE_MOSAIC:
        mosaic_path.unlink()

    run_command(
        [
            "gdal_translate",
            "-of",
            "GTiff",
            "-ot",
            "Byte",
            "-a_nodata",
            OUTPUT_NODATA,
            "-co",
            "TILED=YES",
            "-co",
            "BLOCKXSIZE=512",
            "-co",
            "BLOCKYSIZE=512",
            "-co",
            "COMPRESS=DEFLATE",
            "-co",
            "ZLEVEL=6",
            "-co",
            "PREDICTOR=1",
            "-co",
            "BIGTIFF=YES",
            vrt_path,
            mosaic_path,
        ]
    )
    print(f"Created annual mosaic: {mosaic_path}")


def create_patch_buffer(mosaic_path, final_output_path):
    """Filter vegetation patches, create their buffer, and save the result."""
    if final_output_path.exists() and not OVERWRITE_FINAL_OUTPUT:
        print(f"Final output exists; skipping: {final_output_path}")
        return

    start = perf_counter()
    print(f"Reading annual mosaic: {mosaic_path}")
    with rio.open(mosaic_path) as source:
        profile = source.profile.copy()
        binary_mosaic = source.read(1)
        mosaic_nodata = source.nodata

    valid = (
        np.ones(binary_mosaic.shape, dtype=bool)
        if mosaic_nodata is None
        else binary_mosaic != mosaic_nodata
    )
    high_coverage = binary_mosaic == 1
    del binary_mosaic
    gc.collect()

    print(f"High-coverage pixels before filtering: {high_coverage.sum():,}")

    # Eight-neighbor connectivity includes horizontal, vertical, and diagonal
    # connections between high-coverage pixels.
    structure = generate_binary_structure(rank=2, connectivity=2)
    labeled, patch_count = label(high_coverage, structure=structure)
    del high_coverage
    gc.collect()
    print(f"Connected patches: {patch_count:,}")

    pixel_counts = np.bincount(labeled.ravel())
    retained_labels = np.flatnonzero(pixel_counts >= MIN_PATCH_PIXELS)
    retained_labels = retained_labels[retained_labels != 0]
    retained_patch_count = int(retained_labels.size)
    retained = (
        np.isin(labeled, retained_labels)
        if retained_labels.size
        else np.zeros(labeled.shape, dtype=bool)
    )
    del labeled, pixel_counts, retained_labels
    gc.collect()

    retained_pixels = int(retained.sum())
    retained_area_km2 = retained_pixels * PIXEL_SIZE_M**2 / 1_000_000
    print(
        f"Retained patches: {retained_patch_count:,}\n"
        f"Retained core: {retained_pixels:,} pixels "
        f"({retained_area_km2:,.2f} km²)"
    )

    if retained_pixels == 0:
        buffered = np.zeros(retained.shape, dtype=bool)
    else:
        distance_m = distance_transform_edt(
            ~retained,
            sampling=(PIXEL_SIZE_M, PIXEL_SIZE_M),
        )
        buffered = distance_m <= BUFFER_DISTANCE_M
        del distance_m
        gc.collect()

    if INCLUDE_VEGETATED_CORE:
        output = buffered.astype(np.uint8)
        output_mode = "plusvegetated"
    else:
        output = (buffered & ~retained).astype(np.uint8)
        output_mode = "ringonly"

    del retained, buffered
    gc.collect()

    output[~valid] = OUTPUT_NODATA
    valid_count = int(valid.sum())
    output_count = int(np.count_nonzero(output == 1))
    percentage = output_count / valid_count * 100 if valid_count > 0 else 0.0
    del valid
    gc.collect()

    print(
        f"Output mode: {output_mode}\n"
        f"Pixels with value 1: {output_count:,}\n"
        f"Percent of valid area: {percentage:.2f}%"
    )

    profile.update(
        dtype=rio.uint8,
        count=1,
        nodata=OUTPUT_NODATA,
        compress="deflate",
        zlevel=6,
        predictor=1,
        tiled=True,
        blockxsize=OUTPUT_BLOCK_SIZE,
        blockysize=OUTPUT_BLOCK_SIZE,
        BIGTIFF="YES",
    )
    if final_output_path.exists() and OVERWRITE_FINAL_OUTPUT:
        final_output_path.unlink()
    with rio.open(final_output_path, "w", **profile) as destination:
        destination.write(output, 1)

    del output
    gc.collect()
    print(
        f"Saved {final_output_path}\n"
        f"Patch and buffer time: {(perf_counter() - start) / 60:.1f} minutes"
    )


def process_year(year):
    """Run the complete workflow for one year."""
    print(f"\nProcessing year {year}", flush=True)
    annual_start = perf_counter()

    chunk_paths = sorted(
        INPUT_CHUNK_DIR.glob(INPUT_CHUNK_PATTERN.format(year=year))
    )
    if len(chunk_paths) != EXPECTED_CHUNKS:
        print(
            f"Expected {EXPECTED_CHUNKS} chunks for {year}, but found "
            f"{len(chunk_paths)}; skipping this year."
        )
        return

    inspect_and_validate_chunks(chunk_paths)

    binary_chunk_paths = []
    for input_path in chunk_paths:
        output_path = BINARY_CHUNK_DIR / (
            f"{input_path.stem}_wv75_binary_uint8.tif"
        )
        reclassify_chunk_to_uint8(
            input_path,
            output_path,
            VEGETATION_COVERAGE_THRESHOLD,
        )
        binary_chunk_paths.append(output_path)

    missing = [path for path in binary_chunk_paths if not path.exists()]
    if missing:
        print("Missing binary chunks; skipping annual mosaic:")
        for path in missing:
            print(f"  {path}")
        return

    vrt_path = BINARY_MOSAIC_DIR / f"CONUS_NLCD_{year}_wv75_binary.vrt"
    mosaic_path = BINARY_MOSAIC_DIR / (
        f"CONUS_NLCD_{year}_wv75_binary_mosaic_uint8.tif"
    )
    build_binary_mosaic(binary_chunk_paths, vrt_path, mosaic_path)
    if not mosaic_path.exists():
        print("Annual mosaic was not created; skipping patch analysis.")
        return

    output_mode = "plusvegetated" if INCLUDE_VEGETATED_CORE else "ringonly"
    final_output_path = FINAL_OUTPUT_DIR / (
        f"CONUS_NLCD_{year}_wv75_patch5km2_buffer2400m_{output_mode}.tif"
    )
    create_patch_buffer(mosaic_path, final_output_path)
    print(
        f"Completed {year} in {(perf_counter() - annual_start) / 60:.1f} minutes"
    )


def main():
    """Validate configuration and process all requested years."""
    require_gdal_commands()
    if not INPUT_CHUNK_DIR.exists():
        raise FileNotFoundError(
            f"Input directory does not exist: {INPUT_CHUNK_DIR}. "
            "Replace INPUT_CHUNK_DIR with the directory containing the inputs."
        )
    for directory in (BINARY_CHUNK_DIR, BINARY_MOSAIC_DIR, FINAL_OUTPUT_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    years = years_to_process()
    print("High-vegetation patch and buffer workflow")
    print(f"Input directory: {INPUT_CHUNK_DIR}")
    print(f"Output directory: {OUTPUT_ROOT}")
    print(f"Coverage threshold: {VEGETATION_COVERAGE_THRESHOLD}")
    print(f"Minimum patch area: {MIN_PATCH_AREA_M2 / 1_000_000:.2f} km²")
    print(f"Minimum patch size: {MIN_PATCH_PIXELS:,} pixels")
    print(f"Buffer distance: {BUFFER_DISTANCE_M:,} m")
    print(f"Years: {', '.join(years)}")

    for year in years:
        process_year(year)

    print("All requested years finished.")


if __name__ == "__main__":
    main()
