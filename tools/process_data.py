#!/usr/bin/env python3

"""Script to process output data from running main, density, reblock and analysis."""

import shutil
from pathlib import Path

import numpy as np

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

# Read needed constants from input file
with Path("input.toml").open("rb") as f:
    input = tomllib.load(f)

N = input["system"]["particles"]
T = input["system"]["temperature"]
DIPOLE_MAGNITUDE = np.linalg.norm(input["system"]["dipole_vector"])
DIPOLE_TILT = np.arcsin(input["system"]["dipole_vector"][1] / DIPOLE_MAGNITUDE)
WALL_COEFFICIENT = input["system"]["wall_coefficient"]
ITERATIONS = input["metropolis"]["iterations"]
TRAP = "Pancake" if input["system"]["trapping_frequency_z"] > input["system"]["trapping_frequency_transverse"] else "Cigar"

if __name__ == "__main__":
    path = (
        Path("Outputs")
        / f"{TRAP}/Dipole_str{DIPOLE_MAGNITUDE:.2e}"
        / f"c6{WALL_COEFFICIENT:.2e}"
        / f"Tilt{DIPOLE_TILT:.2e}"
        / f"N{N:.2e}/T{T:.2e}"
        / f"Iterations{ITERATIONS:.2e}"
    )
    path.mkdir(parents=True)
    output_file = Path("configuration.out")
    if output_file.is_file():
        output_file.unlink()
    shutil.copy("input.toml", path / "input.txt")
    for file in [f for f in Path().iterdir() if "out" in f.name or "png" in f.name]:
        file.rename(path / file)
