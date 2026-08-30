from pathlib import Path
import subprocess
import gc

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


# Settings

vtk_dir = Path("path/to/vtk/files") # Folder containing your .vtk files

output_video = vtk_dir / "CC3D_simulation.mp4"

fps = 7

# Use every nth VTK file
stride = 1

output_size = None # None = keep original VTK resolution


# Cell Colors

COLORS = {
    0: (0, 0, 0),         # Medium     → black
    1: (219, 0, 0),       # Tumour     → red
    2: (255, 255, 255),   # CAF        → white
    3: (208, 104, 252),    # myCAF      → purple
    4: (53, 101, 255),    # CD8 T      → blue
}


# FFMPEG Settings

ffmpeg = "ffmpeg"

crf = "16" # high quality

preset = "slow" # slow = better compression/quality, fast = faster encoding


# VTK Files

vtk_files = sorted(vtk_dir.glob("*.vtk"))

if not vtk_files:
    raise RuntimeError(
        f"No VTK files found in {vtk_dir}"
    )

print(f"Found {len(vtk_files)} VTK files")

# Select every nth file
vtk_files = vtk_files[::stride]

print(
    f"Processing {len(vtk_files)} files "
    f"(stride = {stride})"
)

print(f"Output: {output_video}")

# Function: Resize

def resize_image(image, size):

    if size is None:
        return image

    from PIL import Image

    pil_image = Image.fromarray(image)

    pil_image = pil_image.resize(
        size,
        Image.Resampling.NEAREST
    )

    return np.asarray(pil_image)


# Rread the first file separately to obtain image dimensions before starting FFmpeg.

first_file = vtk_files[0]

print()
print("Reading first VTK:")
print(first_file)


reader = vtk.vtkDataSetReader()
reader.SetFileName(str(first_file))
reader.ReadAllScalarsOn()
reader.Update()

mesh = reader.GetOutput()

dimensions = mesh.GetDimensions()

nx = dimensions[0]
ny = dimensions[1]
nz = dimensions[2]

print(
    f"VTK dimensions: "
    f"{nx} × {ny} × {nz}"
)

if nz != 1:
    raise RuntimeError(
        "This script expects a 2D CC3D simulation "
        f"with nz = 1. Found nz = {nz}"
    )


# Check cell type

cell_data = mesh.GetCellData()
point_data = mesh.GetPointData()

cell_array = cell_data.GetArray("CellType")
point_array = point_data.GetArray("CellType")


if cell_array is not None:

    print("Found CellType in CELL data")

    vtk_array = cell_array

elif point_array is not None:

    print("Found CellType in POINT data")

    vtk_array = point_array

else:

    print()
    print("Available cell arrays:")

    for i in range(cell_data.GetNumberOfArrays()):
        print(
            "   ",
            cell_data.GetArrayName(i)
        )

    print()
    print("Available point arrays:")

    for i in range(point_data.GetNumberOfArrays()):
        print(
            "   ",
            point_data.GetArrayName(i)
        )

    raise RuntimeError(
        "Could not find a CellType array."
    )


# Output size
test_width = nx
test_height = ny


# Start FFMPEG

ffmpeg_command = [
    ffmpeg,

    "-y",

    # Input format
    "-f", "rawvideo",

    # Pixel format coming from NumPy
    "-pix_fmt", "rgb24",

    # Dimensions
    "-s",
    f"{test_width}x{test_height}",

    # Frame rate
    "-r", str(fps),

    # Read frames from stdin
    "-i", "-",

    # Video codec
    "-c:v", "libx264",

    # Quality
    "-crf", crf,

    # Encoding speed
    "-preset", preset,

    # Compatibility
    "-pix_fmt", "yuv420p",

    str(output_video)
]


print()
print("Starting FFmpeg...")
print(" ".join(ffmpeg_command))
print()


process = subprocess.Popen(
    ffmpeg_command,
    stdin=subprocess.PIPE
)


# Process VTK files

for index, vtk_file in enumerate(vtk_files):

    print(
        f"[{index + 1}/{len(vtk_files)}] "
        f"{vtk_file.name}",
        flush=True
    )

    # Read VTK

    reader = vtk.vtkDataSetReader()
    reader.SetFileName(str(vtk_file))
    reader.ReadAllScalarsOn()
    reader.Update()

    mesh = reader.GetOutput()


    # Get cell type

    cell_data = mesh.GetCellData()
    point_data = mesh.GetPointData()

    cell_array = cell_data.GetArray("CellType")
    point_array = point_data.GetArray("CellType")

    if cell_array is not None:
        vtk_array = cell_array

    elif point_array is not None:
        vtk_array = point_array

    else:
        process.stdin.close()
        process.wait()

        raise RuntimeError(
            f"CellType not found in {vtk_file}"
        )


    # Convert from VTK to Numpy

    cell_type = vtk_to_numpy(vtk_array)


    # Reshape

    cell_type = cell_type.reshape(
        (ny, nx),
        order="C"
    )


    # Create RGB image

    rgb = np.zeros(
        (ny, nx, 3),
        dtype=np.uint8
    )


    # Apply colours

    actual_types = np.unique(cell_type)

    for cell_type_value, color in COLORS.items():

        mask = cell_type == cell_type_value

        rgb[mask] = color


    # Report unknown types

    unknown_types = set(
        actual_types.tolist()
    ) - set(COLORS.keys())

    if unknown_types:

        print(
            f"   WARNING: unknown CellType values: "
            f"{sorted(unknown_types)}"
        )


    # Orientation

    rgb = np.flipud(image)(rgb)


    # Resize

    rgb = resize_image(
        rgb,
        output_size
    )


    # Send frames to FFMPEG

    process.stdin.write(
        np.ascontiguousarray(rgb).tobytes()
    )


    # Clean memory

    del mesh
    del cell_type
    del rgb
    del reader

    gc.collect()


# Finish FFMPEG

print()
print("Finishing video...")

process.stdin.close()

return_code = process.wait()


if return_code != 0:

    raise RuntimeError(
        f"FFmpeg failed with return code "
        f"{return_code}"
    )


print()
print("============================================")
print("VIDEO COMPLETE")
print("============================================")
print()
print(f"Saved to:")
print(output_video)
