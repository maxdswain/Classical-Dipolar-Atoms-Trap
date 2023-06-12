#!/usr/bin/env python3

import os
import shutil

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

# Read needed constants from input file
with open("input.toml", "rb") as f:
    input = tomllib.load(f)

N = input["system"]["particles"]
T = input["system"]["temperature"]
ITERATIONS = input["metropolis"]["iterations"]
if input["system"]["trapping_frequency_z"] > input["system"]["trapping_frequency_transverse"]:
    TRAP = "pancake"
else:
    TRAP = "cigar"

if __name__ == "__main__":
    folder_name = f"{TRAP}_N={N}_T={T}_{ITERATIONS}_iterations/"
    path = os.path.join("Outputs/", folder_name)
    os.mkdir(path)
    if os.path.isfile("configuration.out"):
        os.remove("configuration.out")
    shutil.copy("input.toml", os.path.join(path, "input.txt"))
    for file in [f for f in os.listdir(".") if "out" in f or "png" in f]:
        os.rename(file, os.path.join(path, file))
