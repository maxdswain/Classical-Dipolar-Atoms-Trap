#!/usr/bin/env python3

import os
import shutil

import numpy as np

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

# Read needed constants from input file
with open("input.toml", "rb") as f:
    input = tomllib.load(f)

N = input["system"]["particles"]
T = input["system"]["temperature"]
DIPOLE_MAGNITUDE = np.linalg.norm(input["system"]["dipole_vector"])
DIPOLE_TILT = np.arcsin(input["system"]["dipole_vector"][1] / DIPOLE_MAGNITUDE)
WALL_COEFFICIENT = input["system"]["wall_coefficient"]
ITERATIONS = input["metropolis"]["iterations"]
if input["system"]["trapping_frequency_z"] > input["system"]["trapping_frequency_transverse"]:
    TRAP = "Pancake"
else:
    TRAP = "Cigar"

if __name__ == "__main__":
    dir = f"{TRAP}/Dipole_str{DIPOLE_MAGNITUDE:.2e}/c6{WALL_COEFFICIENT:.2e}/Tilt{DIPOLE_TILT:.2e}/N{N:.2e}/T{T:.2e}/Iterations{ITERATIONS:.2e}/"
    path = os.path.join("Outputs/", dir)
    os.makedirs(path)
    if os.path.isfile("configuration.out"):
        os.remove("configuration.out")
    shutil.copy("input.toml", os.path.join(path, "input.txt"))
    for file in [f for f in os.listdir(".") if "out" in f or "png" in f]:
        os.rename(file, os.path.join(path, file))
