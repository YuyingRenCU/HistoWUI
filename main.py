#!/usr/bin/env python3
"""Run the HistoWUI raster workflow using the scripts in this directory.

This entry point connects the published processing steps without changing
their classification criteria. Run ``python main.py --help`` for options.
"""

import argparse
import importlib.util
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_YEARS = tuple(range(1985, 2021, 5))
CONTEMPORARY_YEAR = 2020
RADIUS_PIXELS = 19
PIXEL_SIZE_M = 30
BUILDING_TYPE = "ALL"
PARTS = ("part1", "part2", "part3")


def parse_years(value):
    """Parse a comma-separated list of years, preserving its order."""
    try:
        years = tuple(dict.fromkeys(int(item.strip()) for item in value.split(",")))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "years must be comma-separated integers, e.g. 1990,2000,2020"
        ) from exc
    if not years or any(year < 1900 or year > 2100 for year in years):
        raise argparse.ArgumentTypeError("supply at least one valid four-digit year")
    return years


def command_line():
    parser = argparse.ArgumentParser(
        description="Create HistoWUI maps from prepared 30 m building and NLCD rasters."
    )
    parser.add_argument(
        "--building-parts",
        type=Path,
        required=True,
        help="Directory of per-pixel building-count GeoTIFF parts.",
    )
    parser.add_argument(
        "--nlcd-parts",
        type=Path,
        required=True,
        help="Directory of NLCD GeoTIFF parts for neighborhood vegetation coverage.",
    )
    parser.add_argument(
        "--nlcd-reference",
        type=Path,
        required=True,
        help="Directory of annual, full-extent NLCD reference GeoTIFFs.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Directory for intermediate rasters and final HistoWUI maps.",
    )
    parser.add_argument(
        "--mode",
        choices=("historical", "contemporary", "both"),
        default="historical",
        help="Run the historical series, contemporary 2020, or both (default: historical).",
    )
    parser.add_argument(
        "--years",
        type=parse_years,
        default=DEFAULT_YEARS,
        help="Historical years, comma-separated (default: 1985 to 2020 every five years).",
    )
    parser.add_argument(
        "--nlcd-reference-pattern",
        default="CONUS_NLCD_{year}_reextent.tif",
        help="Reference filename pattern containing {year}.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Workers per neighborhood raster (default: 4).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing intermediate and final rasters.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show the planned files and steps without reading inputs or writing outputs.",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    if "{year}" not in args.nlcd_reference_pattern:
        parser.error("--nlcd-reference-pattern must contain {year}")
    return args


def load_script(module_name, filename):
    """Import a sibling script, including filenames containing spaces."""
    script_path = SCRIPT_DIR / filename
    if not script_path.is_file():
        raise FileNotFoundError(f"Required workflow script is missing: {script_path}")
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def selected_years(args):
    """Return years needed for vegetation and each building inventory."""
    historical = args.years if args.mode in ("historical", "both") else ()
    contemporary = (CONTEMPORARY_YEAR,) if args.mode in ("contemporary", "both") else ()
    vegetation = tuple(sorted(set(historical) | set(contemporary)))
    return historical, contemporary, vegetation


def output_paths(args):
    radius_m = RADIUS_PIXELS * PIXEL_SIZE_M
    root = args.output_root.expanduser().resolve()
    return {
        "root": root,
        "building": root / f"building_density_{radius_m}m_conus_split" / BUILDING_TYPE,
        "vegetation": root / f"Vegetation_coverage_{radius_m}m_conus_split",
        "buffer_root": root / "vegetation_buffer",
        "historical": root / "HistoWUI_time_series",
        "contemporary": root / "HistoWUI_contemporary_2020",
    }


def show_plan(args, paths):
    historical, contemporary, vegetation = selected_years(args)
    print(f"Historical years: {', '.join(map(str, historical)) or 'none'}")
    print(f"Contemporary years: {', '.join(map(str, contemporary)) or 'none'}")
    print(f"Vegetation years: {', '.join(map(str, vegetation))}")
    print(f"Radius: {RADIUS_PIXELS} pixels ({RADIUS_PIXELS * PIXEL_SIZE_M} m)")
    print(f"Building parts: {args.building_parts}")
    print(f"NLCD parts: {args.nlcd_parts}")
    print(f"NLCD references: {args.nlcd_reference}")
    print(f"Building outputs: {paths['building']}")
    print(f"Vegetation outputs: {paths['vegetation']}")
    print(f"Vegetation buffers: {paths['buffer_root'] / 'wv75_patch5km2_buffer2400m'}")
    if historical:
        print(f"Historical HistoWUI: {paths['historical']}")
    if contemporary:
        print(f"Contemporary HistoWUI: {paths['contemporary']}")
    print("Order: building counts and vegetation coverage; vegetation patch buffer; WUI classes")


def verify_input_files(args, building, vegetation, historical, contemporary, years):
    """Check all requested source files before starting long raster operations."""
    for label, directory in (
        ("building parts", args.building_parts),
        ("NLCD parts", args.nlcd_parts),
        ("NLCD references", args.nlcd_reference),
    ):
        if not directory.is_dir():
            raise FileNotFoundError(f"{label} directory does not exist: {directory}")

    vegetation.INPUT_DIR = args.nlcd_parts
    veg_files = sorted(args.nlcd_parts.glob("*.tif"))
    veg_by_year = vegetation.group_files_by_year(veg_files)
    for year in years:
        key = str(year)
        if key not in veg_by_year:
            raise FileNotFoundError(f"No NLCD parts found for {year}")
        vegetation.validate_parts(key, veg_by_year[key])
        reference = args.nlcd_reference / args.nlcd_reference_pattern.format(year=year)
        if not reference.is_file():
            raise FileNotFoundError(f"NLCD reference raster is missing: {reference}")

    building.INPUT_DIR = args.building_parts
    for is_contemporary, requested in ((False, historical), (True, contemporary)):
        if not requested:
            continue
        building.USE_CONTEMPORARY = is_contemporary
        files, _ = building.select_input_files()
        by_year = building.group_files_by_year(files)
        for year in requested:
            key = str(year)
            if key not in by_year:
                inventory = "contemporary" if is_contemporary else "historical"
                raise FileNotFoundError(f"No {inventory} building parts found for {year}")
            building.validate_parts(key, by_year[key])
    return veg_by_year


def run_building_year(building, year, is_contemporary, paths, overwrite):
    """Process one year of building parts using the selected inventory."""
    building.USE_CONTEMPORARY = is_contemporary
    files, prefix = building.select_input_files()
    by_year = building.group_files_by_year(files)
    year_files = by_year[str(year)]
    building.validate_parts(str(year), year_files)
    building.OVERWRITE_OUTPUTS = overwrite
    paths["building"].mkdir(parents=True, exist_ok=True)
    radius_m = RADIUS_PIXELS * PIXEL_SIZE_M
    for input_path in year_files:
        part = building.extract_part_id(input_path)
        output_path = paths["building"] / (
            f"{prefix}_{BUILDING_TYPE}_bdgden_{radius_m}m_{year}_{part}.tif"
        )
        building.process_file(input_path, output_path)


def run_vegetation_year(vegetation, year, year_files, paths, overwrite):
    """Process one year of NLCD parts into neighborhood vegetation coverage."""
    vegetation.OVERWRITE_OUTPUTS = overwrite
    paths["vegetation"].mkdir(parents=True, exist_ok=True)
    radius_m = RADIUS_PIXELS * PIXEL_SIZE_M
    for input_path in year_files:
        part = vegetation.extract_part_id(input_path)
        output_path = paths["vegetation"] / (
            f"conus_wv_coverage_{radius_m}m_{year}_{part}.tif"
        )
        vegetation.process_file(input_path, output_path)


def configure_buffer(buffer_script, paths, overwrite):
    """Connect the vegetation output names to the patch buffer step."""
    buffer_script.INPUT_CHUNK_DIR = paths["vegetation"]
    buffer_script.INPUT_CHUNK_PATTERN = (
        f"conus_wv_coverage_{RADIUS_PIXELS * PIXEL_SIZE_M}m_{{year}}_part*.tif"
    )
    buffer_script.OUTPUT_ROOT = paths["buffer_root"]
    buffer_script.BINARY_CHUNK_DIR = paths["buffer_root"] / "binary_wv75_chunks"
    buffer_script.BINARY_MOSAIC_DIR = paths["buffer_root"] / "binary_wv75_annual_mosaics"
    buffer_script.FINAL_OUTPUT_DIR = (
        paths["buffer_root"] / "wv75_patch5km2_buffer2400m"
    )
    buffer_script.INCLUDE_VEGETATED_CORE = True
    buffer_script.OVERWRITE_BINARY_CHUNKS = overwrite
    buffer_script.OVERWRITE_MOSAIC = overwrite
    buffer_script.OVERWRITE_FINAL_OUTPUT = overwrite
    for directory in (
        buffer_script.BINARY_CHUNK_DIR,
        buffer_script.BINARY_MOSAIC_DIR,
        buffer_script.FINAL_OUTPUT_DIR,
    ):
        directory.mkdir(parents=True, exist_ok=True)


def run_buffer_year(buffer_script, year):
    """Create and confirm the vegetation patch and buffer raster."""
    buffer_script.process_year(str(year))
    expected = buffer_script.FINAL_OUTPUT_DIR / (
        f"CONUS_NLCD_{year}_wv75_patch5km2_buffer2400m_plusvegetated.tif"
    )
    if not expected.is_file():
        raise RuntimeError(f"Vegetation buffer was not created for {year}: {expected}")


def configure_wui(wui, args, paths, inventory):
    """Connect the selected building inventory and other prepared inputs."""
    wui.BUILDING_TYPE = BUILDING_TYPE
    wui.BUILDING_DENSITY_DIR = paths["building"]
    wui.BUILDING_DENSITY_PATTERN = f"{inventory}_*{{year}}*{{part}}.tif"
    wui.VEGETATION_COVERAGE_DIR = paths["vegetation"]
    wui.VEGETATION_COVERAGE_PATTERN = (
        f"conus_wv_coverage_{RADIUS_PIXELS * PIXEL_SIZE_M}m_{{year}}_{{part}}.tif"
    )
    wui.VEGETATION_BUFFER_DIR = paths["buffer_root"] / "wv75_patch5km2_buffer2400m"
    wui.NLCD_DIR = args.nlcd_reference
    wui.NLCD_PATTERN = args.nlcd_reference_pattern
    wui.OUTPUT_DIR = paths["historical" if inventory == "conus" else "contemporary"]
    wui.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def run_wui_year(wui, year, inventory, overwrite):
    """Classify one year and give contemporary 2020 a distinct filename."""
    temporary_name = wui.OUTPUT_DIR / f"HistoWUI_{BUILDING_TYPE}_{year}.tif"
    final_name = (
        temporary_name
        if inventory == "conus"
        else wui.OUTPUT_DIR / f"HistoWUI_{BUILDING_TYPE}_{year}_compl.tif"
    )
    if final_name.exists() and not overwrite:
        print(f"WUI output exists; skipping: {final_name}")
        return
    wui.process_year(str(year))
    if inventory != "conus":
        temporary_name.replace(final_name)
    if not final_name.is_file():
        raise RuntimeError(f"WUI raster was not created: {final_name}")


def main():
    args = command_line()
    paths = output_paths(args)
    historical, contemporary, years = selected_years(args)
    show_plan(args, paths)
    if args.dry_run:
        return

    # Load the processing scripts only for a real run so --help and --dry-run
    # work before scientific Python and GDAL are installed.
    building = load_script("Step3_building_den", "Step3_building_den.py")
    vegetation = load_script("Step3_vegetation_cov", "Step3_vegetation_cov.py")
    buffer_script = load_script(
        "Step2_vegetation_buffer",
        "Step2_identify_wv_cluster_and 24buffer_for_interfaceWUI.py",
    )
    wui = load_script("Step3_identify_WUI", "Step3_identify_WUI.py")

    building.RADIUS_PIXELS = RADIUS_PIXELS
    building.PIXEL_SIZE_M = PIXEL_SIZE_M
    building.NUM_WORKERS = args.workers
    vegetation.RADIUS_PIXELS = RADIUS_PIXELS
    vegetation.PIXEL_SIZE_M = PIXEL_SIZE_M
    vegetation.NUM_WORKERS = args.workers
    building.EXPECTED_PARTS = PARTS
    vegetation.EXPECTED_PARTS = PARTS
    buffer_script.EXPECTED_CHUNKS = len(PARTS)
    wui.PARTS = PARTS

    vegetation_files = verify_input_files(
        args, building, vegetation, historical, contemporary, years
    )
    buffer_script.require_gdal_commands()

    for year in historical:
        run_building_year(building, year, False, paths, args.overwrite)
    for year in contemporary:
        run_building_year(building, year, True, paths, args.overwrite)

    for year in years:
        run_vegetation_year(
            vegetation, year, vegetation_files[str(year)], paths, args.overwrite
        )

    configure_buffer(buffer_script, paths, args.overwrite)
    for year in years:
        run_buffer_year(buffer_script, year)

    for year in historical:
        configure_wui(wui, args, paths, "conus")
        run_wui_year(wui, year, "conus", args.overwrite)
    for year in contemporary:
        configure_wui(wui, args, paths, "contemporary")
        run_wui_year(wui, year, "contemporary", args.overwrite)

    print("HistoWUI workflow complete.")


if __name__ == "__main__":
    main()
