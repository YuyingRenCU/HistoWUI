#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate local wildland-vegetation coverage rasters from NLCD inputs.

The script groups spatial NLCD parts by year and calls
``Step1_Function_density_coverage.main`` for each part. Output values range
from 0 to 10,000, where 10,000 represents 100% wildland vegetation within the
configured circular neighborhood.

Input parts processed independently should include an overlap of at least the
neighborhood radius to avoid edge effects at internal part boundaries.
"""

import re
from pathlib import Path

from Step1_Function_density_coverage import (
    main as calculate_neighborhood_values,
)


# -----------------------------------------------------------------------------
# User configuration
# -----------------------------------------------------------------------------
# Replace these paths with the actual input and output locations.
INPUT_DIR = Path("/path/to/nlcd_raster_parts")
OUTPUT_ROOT = Path("/path/to/vegetation_coverage_outputs")

# Radius is expressed in raster pixels. At 30 m resolution:
#   8 pixels  = 240 m, approximately a 250 m neighborhood radius
#   17 pixels = 510 m, approximately a 500 m neighborhood radius
#   19 pixels = 570 m
#   25 pixels = 750 m
RADIUS_PIXELS = 19
PIXEL_SIZE_M = 30
NUM_WORKERS = 4

# Set to None to process every year found in the input directory. Otherwise,
# provide the years to process, for example ("2020",) or ("2020", "2021").
SELECTED_YEARS = ("2020",)

# The workflow expects three spatial parts per year. Change this tuple if the
# public input data use a different partitioning scheme.
EXPECTED_PARTS = ("part1", "part2", "part3")

# Existing outputs are retained unless this option is True.
OVERWRITE_OUTPUTS = False

YEAR_PATTERN = re.compile(r"(?:19|20)\d{2}")
PART_PATTERN = re.compile(r"part\d+", flags=re.IGNORECASE)


def extract_year(path):
    """Extract the first four-digit year from an input filename."""
    match = YEAR_PATTERN.search(path.name)
    if match is None:
        raise ValueError(f"Could not find a four-digit year in {path.name}")
    return match.group(0)


def extract_part_id(path):
    """Extract a part identifier such as ``part1`` from a filename."""
    matches = PART_PATTERN.findall(path.stem)
    if len(matches) != 1:
        raise ValueError(
            f"Expected one part identifier in {path.name}, but found {matches}"
        )
    return matches[0].lower()


def group_files_by_year(files):
    """Group NLCD input paths by their filename year."""
    grouped = {}
    for path in files:
        grouped.setdefault(extract_year(path), []).append(path)
    return {year: sorted(paths) for year, paths in sorted(grouped.items())}


def select_years(files_by_year):
    """Apply the optional year selection and validate requested years."""
    if SELECTED_YEARS is None:
        return files_by_year

    requested = [str(year) for year in SELECTED_YEARS]
    missing = [year for year in requested if year not in files_by_year]
    if missing:
        raise FileNotFoundError(
            f"No NLCD input rasters were found for requested years: {missing}"
        )
    return {year: files_by_year[year] for year in requested}


def validate_parts(year, paths):
    """Confirm that each year contains the expected spatial parts once."""
    part_ids = [extract_part_id(path) for path in paths]
    if len(part_ids) != len(set(part_ids)):
        raise ValueError(f"Duplicate part identifiers found for {year}: {part_ids}")

    missing = sorted(set(EXPECTED_PARTS) - set(part_ids))
    unexpected = sorted(set(part_ids) - set(EXPECTED_PARTS))
    if missing or unexpected:
        raise ValueError(
            f"Unexpected part collection for {year}. "
            f"Missing: {missing or 'none'}; unexpected: {unexpected or 'none'}"
        )


def process_file(input_path, output_path):
    """Calculate neighborhood vegetation coverage for one NLCD raster part."""
    if output_path.exists() and not OVERWRITE_OUTPUTS:
        print(f"Output exists; skipping: {output_path}")
        return

    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")
    calculate_neighborhood_values(
        input_path,
        output_path,
        input_type="nlcd",
        num_workers=NUM_WORKERS,
        radius=RADIUS_PIXELS,
        output_dtype="int16",
        pixel_size_m=PIXEL_SIZE_M,
    )


def main():
    """Process all selected NLCD raster parts."""
    if not INPUT_DIR.exists():
        raise FileNotFoundError(
            f"Input directory does not exist: {INPUT_DIR}. "
            "Replace INPUT_DIR with the directory containing the NLCD inputs."
        )

    input_files = sorted(INPUT_DIR.glob("*.tif"))
    if not input_files:
        raise FileNotFoundError(f"No GeoTIFF inputs were found in {INPUT_DIR}")

    files_by_year = select_years(group_files_by_year(input_files))
    radius_m = RADIUS_PIXELS * PIXEL_SIZE_M
    output_dir = OUTPUT_ROOT / f"Vegetation_coverage_{radius_m}m_conus_split"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Radius: {RADIUS_PIXELS} pixels ({radius_m} m)")
    print(f"Workers: {NUM_WORKERS}")
    print(f"Years: {', '.join(files_by_year)}")

    for year, year_files in files_by_year.items():
        validate_parts(year, year_files)
        for input_path in year_files:
            part_id = extract_part_id(input_path)
            output_path = output_dir / (
                f"conus_wv_coverage_{radius_m}m_{year}_{part_id}.tif"
            )
            process_file(input_path, output_path)

    print("Processing complete.")


if __name__ == "__main__":
    main()
