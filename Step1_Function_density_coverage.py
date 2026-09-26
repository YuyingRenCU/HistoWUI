#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculate local building counts or vegetation coverage for a raster.

The module applies a circular moving window to a single-band raster:

* ``input_type="building"`` sums building counts within the neighborhood.
* ``input_type="nlcd"`` calculates the proportion of selected wildland
  vegetation classes and scales the result from 0 to 10,000.

Processing is performed in tiled windows with a padded read around each output
window. Values outside the raster and input NoData values contribute zero to
the neighborhood calculation, matching the original HistoWUI workflow.

When spatial tiles are processed independently, each input tile should include
an overlap of at least ``radius`` pixels. Otherwise, values near internal tile
boundaries will be calculated as though neighboring tiles contain zeros.
"""

import concurrent.futures
import threading
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import convolve


# NLCD classes treated as wildland vegetation.
WILDLAND_VEGETATION_CLASSES = (41, 42, 43, 51, 52, 71, 90, 95)

# Input and output conventions used by the HistoWUI workflow.
BUILDING_NODATA = -999
OUTPUT_NODATA = -999
VEGETATION_SCALE = 10_000
DEFAULT_PIXEL_SIZE_M = 30
DEFAULT_BLOCK_SIZE = 512

VALID_INPUT_TYPES = {"building", "nlcd"}


def create_convolution_kernel(radius, pixel_size_m=DEFAULT_PIXEL_SIZE_M):
    """Create a circular kernel and report its pixel count and area.

    Parameters
    ----------
    radius : int
        Neighborhood radius in pixels.
    pixel_size_m : float, optional
        Raster pixel size in meters.

    Returns
    -------
    kernel : numpy.ndarray
        Binary circular convolution kernel.
    kernel_pixels : int
        Number of pixels included in the kernel.
    neighborhood_area_km2 : float
        Area represented by the included pixels.
    """
    if not isinstance(radius, (int, np.integer)) or radius < 0:
        raise ValueError("radius must be a non-negative integer")
    if pixel_size_m <= 0:
        raise ValueError("pixel_size_m must be greater than zero")

    rows, columns = np.ogrid[-radius : radius + 1, -radius : radius + 1]
    kernel = ((columns**2 + rows**2) <= radius**2).astype(np.int16)
    kernel_pixels = int(kernel.sum())
    neighborhood_area_km2 = (
        kernel_pixels * pixel_size_m**2 / 1_000_000
    )
    return kernel, kernel_pixels, neighborhood_area_km2


def apply_neighborhood_operation(input_type, block, kernel, kernel_pixels):
    """Apply the selected neighborhood operation to a padded raster block."""
    if input_type == "nlcd":
        vegetation = np.isin(
            block,
            WILDLAND_VEGETATION_CLASSES,
        ).astype(np.int16)
        vegetation_count = convolve(
            vegetation,
            kernel,
            mode="constant",
            cval=0,
        )
        return vegetation_count / kernel_pixels * VEGETATION_SCALE

    # Building inputs contain the number of buildings per pixel. NoData values
    # contribute zero to the neighborhood sum.
    building_count = block.astype(np.int32, copy=True)
    building_count[building_count == BUILDING_NODATA] = 0
    return convolve(
        building_count,
        kernel.astype(np.int32),
        mode="constant",
        cval=0,
    )


def padded_window(window, raster_width, raster_height, radius):
    """Return a clipped read window with convolution padding."""
    column_offset = max(0, int(window.col_off) - radius)
    row_offset = max(0, int(window.row_off) - radius)
    column_end = min(
        raster_width,
        int(window.col_off + window.width) + radius,
    )
    row_end = min(
        raster_height,
        int(window.row_off + window.height) + radius,
    )
    return rasterio.windows.Window(
        column_offset,
        row_offset,
        column_end - column_offset,
        row_end - row_offset,
    )


def main(
    infile,
    outfile,
    input_type,
    num_workers,
    radius=19,
    output_dtype="int16",
    pixel_size_m=DEFAULT_PIXEL_SIZE_M,
):
    """Process a single-band raster using a circular moving window.

    Parameters
    ----------
    infile : str or pathlib.Path
        Input raster. Building rasters must contain per-pixel building counts;
        NLCD rasters must use standard NLCD class values.
    outfile : str or pathlib.Path
        Output GeoTIFF path.
    input_type : {"building", "nlcd"}
        Neighborhood operation to apply.
    num_workers : int
        Number of worker threads.
    radius : int, optional
        Circular neighborhood radius in pixels. With 30 m pixels, the default
        radius of 19 corresponds to 570 m.
    output_dtype : str or numpy dtype, optional
        Output raster data type. ``int16`` preserves the original workflow.
    pixel_size_m : float, optional
        Expected square pixel size in meters.

    Notes
    -----
    NLCD vegetation coverage is truncated when converted to an integer output,
    matching the original implementation. The output NoData metadata is -999,
    but the neighborhood operation produces values for every output pixel.
    """
    if input_type not in VALID_INPUT_TYPES:
        raise ValueError(
            f"input_type must be one of {sorted(VALID_INPUT_TYPES)}, "
            f"not {input_type!r}"
        )
    if not isinstance(num_workers, (int, np.integer)) or num_workers < 1:
        raise ValueError("num_workers must be a positive integer")

    infile = Path(infile)
    outfile = Path(outfile)
    if not infile.exists():
        raise FileNotFoundError(f"Input raster does not exist: {infile}")
    outfile.parent.mkdir(parents=True, exist_ok=True)

    kernel, kernel_pixels, neighborhood_area_km2 = create_convolution_kernel(
        radius,
        pixel_size_m,
    )
    print(f"Input type: {input_type}")
    print(f"Radius: {radius} pixels ({radius * pixel_size_m:g} m)")
    print(f"Kernel pixels: {kernel_pixels:,}")
    print(f"Kernel area: {neighborhood_area_km2:.4f} km²")

    with rasterio.open(infile) as source:
        if source.count != 1:
            raise ValueError(
                f"Expected a single-band raster, but {infile} has "
                f"{source.count} bands"
            )

        source_resolution = (abs(source.transform.a), abs(source.transform.e))
        if not np.allclose(source_resolution, (pixel_size_m, pixel_size_m)):
            raise ValueError(
                f"Expected {pixel_size_m:g} m square pixels, but {infile} has "
                f"resolution {source_resolution}"
            )

        profile = source.profile.copy()
        profile.update(
            driver="GTiff",
            count=1,
            dtype=output_dtype,
            nodata=OUTPUT_NODATA,
            compress="LZW",
            tiled=True,
            blockxsize=DEFAULT_BLOCK_SIZE,
            blockysize=DEFAULT_BLOCK_SIZE,
            BIGTIFF="IF_NEEDED",
        )

        with rasterio.open(outfile, "w", **profile) as destination:
            windows = [window for _, window in destination.block_windows(1)]
            print(
                f"Processing {len(windows):,} blocks with "
                f"{num_workers} worker(s)"
            )

            # Rasterio dataset handles are protected because reads and writes
            # to shared handles are not thread-safe.
            read_lock = threading.Lock()
            write_lock = threading.Lock()

            def process_window(window):
                read_window = padded_window(
                    window,
                    source.width,
                    source.height,
                    radius,
                )
                with read_lock:
                    padded_array = source.read(1, window=read_window)

                convolved = apply_neighborhood_operation(
                    input_type,
                    padded_array,
                    kernel,
                    kernel_pixels,
                )
                start_row = int(window.row_off - read_window.row_off)
                start_column = int(window.col_off - read_window.col_off)
                result = convolved[
                    start_row : start_row + int(window.height),
                    start_column : start_column + int(window.width),
                ].astype(output_dtype)

                with write_lock:
                    destination.write(result, 1, window=window)

            with concurrent.futures.ThreadPoolExecutor(
                max_workers=num_workers
            ) as executor:
                futures = [
                    executor.submit(process_window, window) for window in windows
                ]
                for completed, future in enumerate(
                    concurrent.futures.as_completed(futures),
                    start=1,
                ):
                    # Re-raise worker exceptions so an incomplete output is not
                    # reported as successfully processed.
                    future.result()
                    if completed % 10 == 0 or completed == len(futures):
                        percentage = completed / len(futures) * 100
                        print(
                            f"Completed {percentage:.1f}% "
                            f"({completed:,}/{len(futures):,})"
                        )

    print(f"Saved: {outfile}")


if __name__ == "__main__":
    raise SystemExit(
        "This file defines reusable functions. Import main() from a driver "
        "script and provide the input and output paths."
    )
