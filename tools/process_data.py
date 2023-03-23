#!/usr/bin/env python3

import os
import shutil

if __name__ == "__main__":
    folder_name = "4pi_over8/"
    path = os.path.join("Outputs/N=15_rotations/", folder_name)
    os.mkdir(path)
    if os.path.isfile("configuration.out"):
        os.remove("configuration.out")
    shutil.copy("input.toml", os.path.join(path, "input.txt"))
    for file in [f for f in os.listdir(".") if "out" in f or "png" in f]:
        os.rename(file, os.path.join(path, file))
