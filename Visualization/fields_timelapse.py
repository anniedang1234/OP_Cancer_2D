from pathlib import Path
import subprocess
import gc

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize


# Settings

vtk_dir = Path(
    "path/to/vtk/dir"
)

fps = 7

stride = 1

# Fields

FIELDS = {

    "IFN_gamma": {
        "title": "IFN-γ",
        "output": "IFN_gamma.mp4",
        "cmap": "viridis",
        "scale": "linear",
        "units": "",
    },

    "TGF_beta": {
        "title": "TGF-β",
        "output": "TGF_beta.mp4",
        "cmap": "plasma",
        "scale": "linear",
        "units": "",
    },

    "Collagen": {
        "title": "Collagen",
        "output": "Collagen.mp4",
        "cmap": "viridis",
        "scale": "log",
        "units": "",
    },
}


# Find VTK files

vtk_files = sorted(vtk_dir.glob("*.vtk"))

if not vtk_files:

    raise RuntimeError(
        f"No VTK files found in {vtk_dir}"
    )

print(
    f"Found {len(vtk_files)} VTK files"
)

vtk_files = vtk_files[::stride]

print(
    f"Using {len(vtk_files)} files "
    f"(stride={stride})"
)


# Read field

def read_field(filename, field_name):

    reader = vtk.vtkDataSetReader()

    reader.SetFileName(
        str(filename)
    )

    reader.ReadAllScalarsOn()
    reader.Update()

    mesh = reader.GetOutput()

    dimensions = mesh.GetDimensions()

    nx = dimensions[0]
    ny = dimensions[1]
    nz = dimensions[2]

    if nz != 1:

        raise RuntimeError(
            f"{filename} is not 2D: "
            f"{dimensions}"
        )

    point_data = mesh.GetPointData()

    array = point_data.GetArray(
        field_name
    )

    if array is None:

        raise RuntimeError(
            f"Could not find "
            f"{field_name} in {filename}"
        )

    values = vtk_to_numpy(array)

    values = values.reshape(
        (ny, nx),
        order="C"
    )

    values = np.flipud(values)

    return values, nx, ny


# Pass through once to determine global mins & maxes

print()
print("============================================")
print("PASS 1: FINDING GLOBAL FIELD RANGES")
print("============================================")
print()


global_ranges = {}


for field_name in FIELDS:

    print(
        f"Scanning {field_name}..."
    )

    global_min = np.inf
    global_max = -np.inf

    for i, vtk_file in enumerate(vtk_files):

        values, nx, ny = read_field(
            vtk_file,
            field_name
        )

        finite_values = values[
            np.isfinite(values)
        ]

        if finite_values.size == 0:
            continue

        current_min = finite_values.min()
        current_max = finite_values.max()

        global_min = min(
            global_min,
            current_min
        )

        global_max = max(
            global_max,
            current_max
        )

        del values
        gc.collect()

    global_ranges[field_name] = (
        global_min,
        global_max
    )

    print(
        f"   min = {global_min:.8e}"
    )

    print(
        f"   max = {global_max:.8e}"
    )

    print()


# Print final ranges

print()
print("============================================")
print("GLOBAL RANGES")
print("============================================")

for field_name, (
    minimum,
    maximum
) in global_ranges.items():

    print(
        f"{field_name:12s}: "
        f"{minimum:.8e} → {maximum:.8e}"
    )

print()


# Create each video

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

    # Log scale

    if settings["scale"] == "log":

        # LogNorm cannot use zero.
        positive_values = minimum

        if positive_values <= 0:

            # Find smallest positive value
            positive_values = np.finfo(
                np.float32
            ).tiny

        norm = LogNorm(
            vmin=positive_values,
            vmax=maximum
        )

    # Linear scale

    else:

        norm = Normalize(
            vmin=minimum,
            vmax=maximum
        )


    # Output file

    output_file = (
        vtk_dir / settings["output"]
    )


    # Figure settings

    fig, ax = plt.subplots(
        figsize=(9.70, 16.90),
        dpi=100
    )

    fig.subplots_adjust(
        left=0.03,
        right=0.84,
        top=0.94,
        bottom=0.03
    )

    ax.set_axis_off()

    ax.set_title(
        settings["title"],
        fontsize=22,
        pad=12
    )


    # First frame

    values, nx, ny = read_field(
        vtk_files[0],
        field_name
    )

    image = ax.imshow(
        values,
        cmap=settings["cmap"],
        norm=norm,
        origin="upper",
        interpolation="nearest"
    )


    # Colour bar

    colorbar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.046,
        pad=0.04
    )

    colorbar.set_label(
        settings["title"],
        fontsize=16
    )

    colorbar.ax.tick_params(
        labelsize=11
    )


    # Time label

    time_text = ax.text(
        0.02,
        0.98,
        vtk_files[0].stem,
        transform=ax.transAxes,
        color="white",
        fontsize=14,
        verticalalignment="top",
        bbox=dict(
            facecolor="black",
            alpha=0.6,
            pad=4
        )
    )


    # Figure size

    fig.canvas.draw()

    width, height = (
        fig.canvas.get_width_height()
    )

    print(
        f"Video resolution: "
        f"{width} × {height}"
    )


    # Start FFMPEG

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
        f"968x1688",

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

    for i, vtk_file in enumerate(
        vtk_files
    ):

        print(
            f"[{i + 1}/{len(vtk_files)}] "
            f"{vtk_file.name}",
            flush=True
        )

        values, nx, ny = read_field(
            vtk_file,
            field_name
        )

        # Update image
        image.set_data(values)

        # Update timestep label
        time_text.set_text(
            vtk_file.stem
        )

        # Render figure
        fig.canvas.draw()

        # Get RGB pixels
        frame = np.asarray(
            fig.canvas.buffer_rgba()
        )[:, :, :3]
        
        frame = frame[:1688, :968, :]

        # Ensure contiguous memory
        frame = np.ascontiguousarray(
            frame,
            dtype=np.uint8
        )

        # Send to FFmpeg
        process.stdin.write(
            frame.tobytes()
        )

        del values
        del frame

        gc.collect()


    # Finish FFMPEG

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
print("ALL FIELD VIDEOS COMPLETE")
print("============================================")
print()

for field_name, settings in FIELDS.items():

    print(
        vtk_dir / settings["output"]
    )
