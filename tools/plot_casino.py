#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np

stop = False
with open("expval.data", "r") as f:
    for line in f:
        if stop:
            grid_size = [int(x) for x in line.rsplit()]
            break
        if line == "Grid size\n":
            stop = True
BINS_X, BINS_Y = grid_size[:2]


def plot_density_contour() -> None:
    plt.rcParams.update({"font.size": 24})
    fig, ax = plt.subplots(figsize=(15, 12))
    if not os.path.isfile("2Dplot.dat"):
        os.system("plot_expval")
    density = np.delete(np.loadtxt("2Dplot.dat"), 2, 1).reshape(BINS_X, BINS_Y, 3)
    cs = ax.contourf(
        density[:, :, 0],
        density[:, :, 1],
        density[:, :, 2],
        40,
        cmap="inferno",
    )
    ax.set(xlabel="$x$ position", ylabel="$y$ position")
    fig.colorbar(cs).ax.set_ylabel("Number density $(a.u.)$", rotation=270, labelpad=30)
    plt.savefig("density_contour.png")


if __name__ == "__main__":
    plot_density_contour()
