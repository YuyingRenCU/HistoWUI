#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate local building-count rasters from per-pixel building inputs.

The script selects either the contemporary building inventory or the historical
time-series inputs, groups files by year, and calls ``neighborhood_convolution``
for each spatial part. Output values are the number of buildings within the
configured circular neighborhood.

Input parts processed independently should include an overlap of at least the
neighborhood radius. For example, a radius of 19 pixels requires a 19-pixel
overlap to prevent edge effects at internal part boundaries.
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
INPUT_DIR = Path("/path/to/building_count_per_pixel_parts")
OUTPUT_ROOT = Path("/path/to/building_neighborhood_outputs")

BUILDING_TYPE = "ALL"
USE_CONTEMPORARY = True

# Radius is expressed in raster pixels. At 30 m resolution:
#   8 pixels  = 240 m, approximately a 250 m neighborhood radius
#   17 pixels = 510 m, approximately a 500 m neighborhood radius
#   19 pixels = 570 m
#   25 pixels = 750 m
RADIUS_PIXELS = 19
PIXEL_SIZE_M = 30
NUM_WORKERS = 4

# The workflow expects three spatial parts per year. Change this tuple if the
# public input data use a different partitioning scheme.
EXPECTED_PARTS = ("part1", "part2", "part3")

# Existing outputs are retained unless this option is True.
OVERWRITE_OUTPUTS = False

YEAR_PATTERN = re.compile(r"(?:19|20)\d{2}")


def extract_year(path):
    """Extract the first four-digit year from an input filename."""
    match = YEAR_PATTERN.search(path.name)
    if match is None:
        raise ValueError(f"Could not find a four-digit year in {path.name}")
    return match.group(0)


def extract_part_id(path):
    """Extract and validate the final underscore-delimited part identifier."""
    part_id = path.stem.split("_")[-1]
    if not re.fullmatch(r"part\d+", part_id, flags=re.IGNORECASE):
        raise ValueError(
            f"Could not identify a part suffix such as 'part1' in {path.name}"
        )
    return part_id.lower()


def select_input_files():
    """Select contemporary or historical input rasters."""
    if USE_CONTEMPORARY:
        files = sorted(INPUT_DIR.glob("contemporary_*.tif"))
        output_prefix = "contemporary"
    else:
        files = sorted(
            path
            for path in INPUT_DIR.glob("*.tif")
            if not path.name.startswith("contemporary_")
        )
        output_prefix = "conus"

    if not files:
        mode = "contemporary" if USE_CONTEMPORARY else "historical"
        raise FileNotFoundError(
            f"No {mode} GeoTIFFs were found in {INPUT_DIR}. "
            "Update INPUT_DIR or the filename pattern."
        )
    return files, output_prefix


def group_files_by_year(files):
    """Group input paths by the year encoded in their filenames."""
    grouped = {}
    for path in files:
        grouped.setdefault(extract_year(path), []).append(path)
    return {year: sorted(paths) for year, paths in sorted(grouped.items())}


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
    """Calculate neighborhood building counts for one raster part."""
    if output_path.exists() and not OVERWRITE_OUTPUTS:
        print(f"Output exists; skipping: {output_path}")
        return

    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")
    calculate_neighborhood_values(
        input_path,
        output_path,
        input_type="building",
        num_workers=NUM_WORKERS,
        radius=RADIUS_PIXELS,
        output_dtype="int16",
        pixel_size_m=PIXEL_SIZE_M,
    )


def main():
    """Process all selected building-count rasters."""
    if not INPUT_DIR.exists():
        raise FileNotFoundError(
            f"Input directory does not exist: {INPUT_DIR}. "
            "Replace INPUT_DIR with the directory containing the inputs."
        )

    radius_m = RADIUS_PIXELS * PIXEL_SIZE_M
    # Retain the HistoWUI filename convention used by downstream scripts.
    # The raster values are neighborhood building totals; "bdgden" is the
    # legacy project label for these locally aggregated values.
    output_dir = (
        OUTPUT_ROOT
        / f"building_density_{radius_m}m_conus_split"
        / BUILDING_TYPE
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    files, output_prefix = select_input_files()
    files_by_year = group_files_by_year(files)

    print(f"Mode: {'contemporary' if USE_CONTEMPORARY else 'historical'}")
    print(f"Radius: {RADIUS_PIXELS} pixels ({radius_m} m)")
    print(f"Files selected: {len(files):,}")
    print(f"Years: {', '.join(files_by_year)}")

    for year, year_files in files_by_year.items():
        validate_parts(year, year_files)
        for input_path in year_files:
            part_id = extract_part_id(input_path)
            output_path = output_dir / (
                f"{output_prefix}_{BUILDING_TYPE}_bdgden_"
                f"{radius_m}m_{year}_{part_id}.tif"
            )
            process_file(input_path, output_path)

    print("Processing complete.")


if __name__ == "__main__":
    main()
