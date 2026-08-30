from pathlib import Path
import subprocess
import numpy as np
import vtk

from vtk.util.numpy_support import vtk_to_numpy

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm


# ============================================================
# SETTINGS
# ============================================================

vtk_dir = Path(
    "path/to/vtk/dir"
)

fps = 7

stride = 1


# Cell colours

CELL_COLORS = {

    0: (0, 0, 0),         # Medium     → not visualized
    1: (219, 0, 0),       # Tumour     → red
    2: (255, 255, 255),   # CAF        → white
    3: (180, 74, 255),    # myCAF      → purple
    4: (53, 101, 255),    # CD8 T      → blue

}


# Fields

FIELDS = {

    "IFN_gamma": {

        "title": "IFN-γ",

        "output": "Cells_IFNgamma.mp4",

        "cmap": "bone",

        "scale": "linear",

        "label": "IFN-γ concentration",

    },

    "TGF_beta": {

        "title": "TGF-β",

        "output": "Cells_TGFbeta.mp4",

        "cmap": "bone",

        "scale": "linear",

        "label": "TGF-β concentration",

    },

    "Collagen": {

        "title": "Collagen",

        "output": "Cells_Collagen.mp4",

        "cmap": "bone",

        "scale": "log",

        "label": "Collagen concentration",

    },

}


# Find VTK files

vtk_files = sorted(
    vtk_dir.glob("Step_*.vtk")
)

if not vtk_files:

    raise RuntimeError(
        f"No Step_*.vtk files found in {vtk_dir}"
    )

vtk_files = vtk_files[::stride]

print(
    f"Found {len(vtk_files)} VTK files"
)

print(
    f"Using every {stride} file"
)


# Read cells & fields

def read_data(filename, field_name):

    reader = vtk.vtkDataSetReader()

    reader.SetFileName(
        str(filename)
    )

    reader.ReadAllScalarsOn()
    reader.Update()

    mesh = reader.GetOutput()

    nx, ny, nz = mesh.GetDimensions()

    if nz != 1:

        raise RuntimeError(
            f"{filename} is not 2D: "
            f"{nx} x {ny} x {nz}"
        )

    point_data = mesh.GetPointData()

    cell_array = point_data.GetArray(
        "CellType"
    )

    field_array = point_data.GetArray(
        field_name
    )

    if cell_array is None:

        raise RuntimeError(
            f"CellType not found in "
            f"{filename}"
        )

    if field_array is None:

        raise RuntimeError(
            f"{field_name} not found in "
            f"{filename}"
        )

    cells = vtk_to_numpy(
        cell_array
    )

    field = vtk_to_numpy(
        field_array
    )

    cells = cells.reshape(
        (ny, nx),
        order="C"
    )

    field = field.reshape(
        (ny, nx),
        order="C"
    )

    # Flip image upside down to match orientation of original data

    cells = np.flipud(cells)

    field = np.flipud(field)

    return cells, field, nx, ny


# Find global field range

print()
print("============================================")
print("FINDING GLOBAL FIELD RANGES")
print("============================================")
print()


global_ranges = {}


for field_name in FIELDS:

    print(
        f"Scanning {field_name}..."
    )

    minimum = np.inf
    maximum = -np.inf

    for filename in vtk_files:

        cells, field, nx, ny = read_data(
            filename,
            field_name
        )

        valid = field[
            np.isfinite(field)
        ]

        if valid.size == 0:
            continue

        minimum = min(
            minimum,
            valid.min()
        )

        maximum = max(
            maximum,
            valid.max()
        )

    global_ranges[field_name] = (
        minimum,
        maximum
    )

    print(
        f"  Minimum: {minimum:.8e}"
    )

    print(
        f"  Maximum: {maximum:.8e}"
    )

    print()


# Create videos

for field_name, settings in FIELDS.items():

    print()
    print("============================================")
    print(
        f"CREATING {settings['title']} VIDEO"
    )
    print("============================================")
    print()

    minimum, maximum = (
        global_ranges[field_name]
    )


    # Normalize colours

    if settings["scale"] == "log":

        # LogNorm cannot use zero.
        #
        # Find the smallest positive value.
        # ----------------------------------------------------

        positive_min = np.inf

        for filename in vtk_files:

            cells, field, nx, ny = read_data(
                filename,
                field_name
            )

            positive = field[
                field > 0
            ]

            if positive.size:

                positive_min = min(
                    positive_min,
                    positive.min()
                )

        norm = LogNorm(
            vmin=positive_min,
            vmax=maximum
        )

        print(
            f"Log scale:"
        )

        print(
            f"  vmin = {positive_min:.8e}"
        )

        print(
            f"  vmax = {maximum:.8e}"
        )

    else:

        norm = Normalize(
            vmin=minimum,
            vmax=maximum
        )

        print(
            f"Linear scale:"
        )

        print(
            f"  vmin = {minimum:.8e}"
        )

        print(
            f"  vmax = {maximum:.8e}"
        )


    # Color map

    cmap = plt.get_cmap(
        settings["cmap"]
    )


    # First frame

    cells, field, nx, ny = read_data(
        vtk_files[0],
        field_name
    )


    # Create background

    background = (
        cmap(norm(field))[:, :, :3] * 255
    ).astype(np.uint8)


    # Layer cells on top of background

    rgb = background.copy()

    for cell_type, color in CELL_COLORS.items():

        if cell_type == 0: # Skip medium cell type
            continue

        mask = (
            cells == cell_type
        )

        rgb[mask] = color


    # Figure

    fig, ax = plt.subplots(
        figsize=(9.70, 16.90),
        dpi=100
    )

    fig.subplots_adjust(
        left=0.02,
        right=0.84,
        top=0.94,
        bottom=0.02
    )


    image = ax.imshow(
        rgb,
        origin="upper",
        interpolation="nearest"
    )

    ax.set_axis_off()


    # Title

    ax.set_title(
        settings["title"] + " + Cells",
        fontsize=20,
        pad=10
    )


    # Colour bar

    scalar_map = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=norm
    )

    scalar_map.set_array([])

    colorbar = fig.colorbar(
        scalar_map,
        ax=ax,
        fraction=0.046,
        pad=0.04
    )

    colorbar.set_label(
        settings["label"],
        fontsize=14
    )


    # Frame / timestep label

    time_text = ax.text(
        0.02,
        0.98,
        vtk_files[0].stem,

        transform=ax.transAxes,

        color="white",

        fontsize=12,

        verticalalignment="top",

        bbox=dict(
            facecolor="black",
            alpha=0.6,
            pad=4
        )
    )


    # Determine output size

    fig.canvas.draw()

    frame = np.asarray(
        fig.canvas.buffer_rgba()
    )[:, :, :3]


    # Make dimensions even

    height = frame.shape[0] - (
        frame.shape[0] % 2
    )

    width = frame.shape[1] - (
        frame.shape[1] % 2
    )

    print(
        f"Video resolution: "
        f"{width} x {height}"
    )


    # FFMPEG

    output_file = (
        vtk_dir / settings["output"]
    )

    ffmpeg_command = [

        "ffmpeg",

        "-y",

        "-f",
        "rawvideo",

        "-vcodec",
        "rawvideo",

        "-pix_fmt",
        "rgb24",

        "-s",
        f"{width}x{height}",

        "-r",
        str(fps),

        "-i",
        "-",

        "-an",

        "-c:v",
        "libx264",

        "-crf",
        "16",

        "-preset",
        "slow",

        "-pix_fmt",
        "yuv420p",

        str(output_file)
    ]


    process = subprocess.Popen(
        ffmpeg_command,
        stdin=subprocess.PIPE
    )


    # Process frames

    for i, filename in enumerate(
        vtk_files
    ):

        print(
            f"[{i + 1}/{len(vtk_files)}] "
            f"{filename.name}",
            flush=True
        )

        cells, field, nx, ny = read_data(
            filename,
            field_name
        )


        # Set background as field

        background = (
            cmap(norm(field))[:, :, :3]
            * 255
        ).astype(np.uint8)


        rgb = background.copy()


        # Set foreground as cells

        for cell_type, color in CELL_COLORS.items():

            if cell_type == 0:
                continue

            mask = (
                cells == cell_type
            )

            rgb[mask] = color


        # Update image

        image.set_data(rgb)

        time_text.set_text(
            filename.stem
        )

        fig.canvas.draw()

        frame = np.asarray(
            fig.canvas.buffer_rgba()
        )[:, :, :3]

        # Crop
        frame = frame[
            :height,
            :width,
            :
        ]

        frame = np.ascontiguousarray(
            frame,
            dtype=np.uint8
        )


        # Send to FFMPEG

        process.stdin.write(
            frame.tobytes()
        )


    # Finish video

    process.stdin.close()

    return_code = process.wait()

    plt.close(fig)


    if return_code != 0:

        raise RuntimeError(
            f"FFmpeg failed for "
            f"{field_name}"
        )


    print()
    print(
        f"Finished: {output_file}"
    )


# Confirm done

print()
print("============================================")
print("ALL THREE VIDEOS COMPLETE")
print("============================================")
print()

print(
    vtk_dir / "Cells_IFNgamma.mp4"
)

print(
    vtk_dir / "Cells_TGFbeta.mp4"
)

print(
    vtk_dir / "Cells_Collagen.mp4"
)
