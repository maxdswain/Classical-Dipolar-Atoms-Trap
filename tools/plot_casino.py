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
plt.rcParams.update({"font.size": 24, "font.family": "Latin Modern Roman", "text.usetex": True, "figure.dpi": 150})


def plot_density_contour() -> None:
    plt.rcParams.update({"font.size": 40})
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
    ax.set(xlabel="$x$ position (a.u.)", ylabel="$y$ position (a.u.)")
    fig.colorbar(cs).ax.set_ylabel("Number density (a.u.)", rotation=270, labelpad=50)
    plt.savefig("density_contour.png")


def plot_densities() -> None:
    plt.rcParams.update({"font.size": 24})
    _, ax = plt.subplots(figsize=(12, 9))
    dipole_moments, qmc = [3.5, 30, 180, 71000], ["VMC", "DMC"]
    for i, density in enumerate(
        map(
            lambda f: np.delete(np.loadtxt(f), 2, 1).reshape(BINS_X, BINS_Y, 3),
            sorted([f for f in os.listdir(".") if os.path.isfile(f) and "2Dplot" in f]),
        )
    ):
        ax.plot(
            density[0, :, 0],
            np.sum(density[:, :, 2], axis=0) / BINS_X,
            label=f"{qmc[0 if i < 4 else 1]} $\ell_{{dd}}={dipole_moments[i % 4]}$",
        )
    ax.set(xlabel="$x$ position (a.u.)", ylabel="Number density (a.u.)")
    ax.legend(loc="upper right")
    plt.savefig("densities.png")


if __name__ == "__main__":
    plot_density_contour()
