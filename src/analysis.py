#!/usr/bin/env python3

"""Functions to read in, process, analyse and then produce plots of the output data from runs of the Classical Metropolis Algorithm."""

from collections.abc import Callable
from itertools import combinations
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import seaborn as sns
from scipy.spatial import distance

# Read needed constants from input file
with Path("input.toml").open("rb") as f:
    input = tomllib.load(f)

N = input["system"]["particles"]
ITERATIONS = input["metropolis"]["iterations"]
T = input["system"]["temperature"]
SAMPLING_RATE = input["metropolis"]["sampling_rate"]
WALL_ORDER = input["system"]["wall_order"]
WALL_COEFFICIENT = input["system"]["wall_coefficient"]
DIPOLE_MOMENT = np.linalg.norm(input["system"]["dipole_vector"])
BINS_X = input["calculations"]["bins_x"]
BINS_Y = input["calculations"]["bins_y"]
BINS_Z = input["calculations"]["bins_z"]
CUTOFF = input["calculations"]["cutoff"]

plt.rcParams.update({"font.size": 34, "font.family": "Latin Modern Roman", "text.usetex": True, "figure.dpi": 150})
round_to_n = lambda x, n: x if x == 0 else round(x, -int(np.floor(np.log10(np.abs(x)))) + (n - 1))


def read_simulation_data() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Read Classical Metropolis Algorithm output data into NumPy arrays."""
    positions = np.loadtxt("position_data.out").reshape(ITERATIONS // SAMPLING_RATE, N, 3)
    energies = np.loadtxt("energy_data.out")
    return positions, energies


def boltzmann_distribution(energies: npt.NDArray[np.float64]) -> Callable[[float], float]:
    """Create Boltzmann distribution function."""
    BOLTZMANN_CONSTANT = 1  # Boltzmann constant in defined systems of units based on values in input.toml
    Z = sum([np.exp(-energy / (BOLTZMANN_CONSTANT * T)) for energy in energies])
    return lambda energy: len(energies) * np.exp(-energy / (BOLTZMANN_CONSTANT * T)) / Z


def plot_energy_histogram(energies: npt.NDArray[np.float64]) -> None:
    """Plot Boltzmann energy histogram."""
    _, ax = plt.subplots()
    sorted_energies = sorted(energies)
    distribution = boltzmann_distribution(sorted_energies)
    ax.hist(energies, bins=20, alpha=0.5)
    ax.plot(
        sorted_energies,
        [distribution(energy) for energy in sorted_energies],
        "r-",
        lw=2,
        label="Boltzmann distribution",
    )
    ax.set(xlabel="Energy", ylabel="Frequency")
    ax.legend(loc="upper right")
    plt.show()


def plot_positions_iterations(positions: npt.NDArray[np.float64], component: int) -> None:
    """Plot positions (y) of the atoms against iterations (x)."""
    _, ax = plt.subplots()
    for i in range(N):
        ax.plot(np.linspace(SAMPLING_RATE, ITERATIONS, ITERATIONS // SAMPLING_RATE), positions[:, i, component])
    ax.set(xlabel="Iterations", ylabel=f"${['x', 'y', 'z'][component]}$ Positions")
    plt.show()


def plot_energies_iterations(energies: npt.NDArray[np.float64]) -> None:
    """Plot total energy (y) of the atoms against iterations (x)."""
    _, ax = plt.subplots(figsize=(15, 9))
    ax.plot(np.linspace(0, ITERATIONS, len(energies)), energies)
    ax.set(xlabel="Iterations", ylabel="Energy", ylim=(0, 5 * round_to_n(energies[-1], 1)))
    plt.savefig("energies_iterations.png")


def plot_many_energies_iterations() -> None:
    """Plot total energy (y) of various runs of the simulation against iterations (x)."""
    _, ax = plt.subplots(figsize=(15, 9))
    energy_data = [np.loadtxt(f) for f in Path().iterdir() if f.is_file() and "energy" in f.name]
    for energies in energy_data:
        ax.plot(np.linspace(0, ITERATIONS, len(energies)), energies)
    ax.set(xlabel="Iterations", ylabel="Energy (a.u.)", xlim=(0, ITERATIONS))
    plt.ylim(1.117 * round_to_n(energy_data[0][-1], 1), 1.124 * round_to_n(energy_data[0][-1], 1))
    plt.savefig("many_energies_iterations.png")


def plot_snapshots(positions: npt.NDArray[np.float64], iteration: int) -> None:
    """Plot a 2D slice of the positions of the atoms at a specified iteration."""
    _, axs = plt.subplots(1, 3, figsize=(15, 4))
    for i, (coord, x) in enumerate(zip(combinations("xyz", 2), combinations(range(3), 2))):
        axs[i].scatter(positions[iteration, :, x[0]], positions[iteration, :, x[1]], c="black")
        axs[i].set(xlabel=f"${coord[0]}$ position", ylabel=f"${coord[1]}$ position")
    plt.tight_layout()
    plt.savefig("snapshots.png")


def plot_error() -> None:
    """Plot error (y) against reblocking transformation number (x)."""
    _, ax = plt.subplots(figsize=(12, 9))
    errors = np.loadtxt("error_data.out")
    ax.plot(np.arange(1, len(errors) + 1), errors, color="black")
    ax.set(xlabel="Reblocking transformation number", ylabel="Standard error", xticks=np.arange(1, len(errors) + 1))
    plt.savefig("error_rtn.png")


def plot_potential() -> None:
    """Plot interatomic potential against distance between two atoms."""
    _, ax = plt.subplots(figsize=(12, 9))
    r = np.linspace(0.9, 3, 101)  # Change ranges of values to ones relevant to system length scales
    repulsive_wall = WALL_COEFFICIENT * r**-WALL_ORDER
    dd_interaction = -2 * DIPOLE_MOMENT**2 * r**-3
    ax.plot(r, dd_interaction + repulsive_wall, label="Combined potential")
    ax.plot(r, repulsive_wall, label=rf"$\frac{{c_{6}}}{{r^{{{WALL_ORDER}}}}}$ potential")
    ax.plot(r, dd_interaction, label="Dipole-dipole interaction")
    ax.set(xlabel="Difference in position, $r$ (a.u.)", ylabel="Potential, $U(r)$ (a.u.)")
    ax.legend(loc="upper right")
    plt.savefig("potential_components.png")


def plot_potentials() -> None:
    """Plot numerous interatomic potentials that vary by repulsion strength against distance between two atoms."""
    _, ax = plt.subplots(figsize=(12, 9))
    r = np.linspace(0.9, 3, 101)  # Change ranges of values to ones relevant to system length scales
    length_scale = [25, 250, 2500, 25000, 250000]
    for i, x in enumerate([WALL_COEFFICIENT * 10**i for i in range(-2, 3)]):
        potential = x * r**-WALL_ORDER - 2 * DIPOLE_MOMENT**2 * r**-3
        ax.plot(r, potential, label=rf"$\ell_{{c_{{6}}}}={length_scale[i]}$")
    ax.set(xlabel="Difference in position, $r$ (a.u.)", ylabel="Potential, $U(r)$ (a.u.)", yscale="symlog")
    ax.legend(loc="upper right")
    plt.savefig("potentials.png")


def plot_density(positions: npt.NDArray[np.float64]) -> None:
    """Plot a colourmap of the number density for a pancake-shaped trap."""
    plt.rcParams.update({"font.size": 40, "xtick.major.pad": 10})
    fig, ax = plt.subplots(figsize=(15, 12))
    density = np.sum(
        [
            np.loadtxt(f).reshape(BINS_X, BINS_Y)
            for f in Path().iterdir()
            if f.is_file() and "density" in f.name and "out" in f.name and "pair" not in f.name
        ],
        axis=0,
    )
    density /= ITERATIONS / SAMPLING_RATE - CUTOFF
    x_bin_length = 1.1 * (np.max(positions[CUTOFF:, :, 0]) - np.min(positions[CUTOFF:, :, 0])) / BINS_X
    y_bin_length = 1.1 * (np.max(positions[CUTOFF:, :, 1]) - np.min(positions[CUTOFF:, :, 1])) / BINS_Y
    cs = ax.contourf(density, 40, cmap="inferno")
    ax.set(
        xlabel="$x$ position (a.u.)",
        ylabel="$y$ position (a.u.)",
        xticks=np.linspace(0, BINS_X - 5, num=7),
        yticks=np.linspace(0, BINS_Y - 5, num=7),
        xticklabels=[round_to_n(x * x_bin_length * BINS_X / 7, 3) for x in range(-3, 4)],
        yticklabels=[round_to_n(y * y_bin_length * BINS_Y / 7, 3) for y in range(-3, 4)],
    )
    fig.colorbar(cs).ax.set_ylabel("Number density (a.u.)", rotation=270, labelpad=50)
    plt.savefig("density_contour.png")


def plot_pair_density(positions: npt.NDArray[np.float64]) -> None:
    """Plot a pair density colourmap for a cigar-shaped trap."""
    plt.rcParams.update({"font.size": 24})
    fig, ax = plt.subplots(figsize=(15, 12))
    pair_density = np.sum(
        [
            np.loadtxt(f).reshape(BINS_Z, BINS_Z)
            for f in Path().iterdir()
            if f.is_file() and "pair" in f.name and "out" in f.name
        ],
        axis=0,
    )
    pair_density /= ITERATIONS / SAMPLING_RATE - CUTOFF
    z_bin_length = 1.1 * (np.max(positions[CUTOFF:, :, 2]) - np.min(positions[CUTOFF:, :, 2])) / BINS_Z
    cs = ax.contourf(pair_density, 40, cmap="inferno")
    ax.set(
        xlabel="$z_{1}$ position (a.u.)",
        ylabel="$z_{2}$ position (a.u.)",
        xticks=np.linspace(0, BINS_Z - 1, num=7),
        yticks=np.linspace(0, BINS_Z - 1, num=7),
        xticklabels=[round_to_n(z * z_bin_length * BINS_Z / 7, 4) for z in range(-3, 4)],
        yticklabels=[round_to_n(z * z_bin_length * BINS_Z / 7, 4) for z in range(-3, 4)],
    )
    fig.colorbar(cs).ax.set_ylabel("Mean number of pairs", rotation=270, labelpad=30)
    plt.savefig("pair_density_contour.png")


def plot_interparticle_distance() -> None:
    """Plot the pairwise interparticle distances between the atoms."""
    _, ax = plt.subplots(figsize=(12, 9))
    pw_size, pos_size = int(0.5 * N * (N - 1)), int(ITERATIONS / SAMPLING_RATE - CUTOFF)
    temps = [1e0, 1e1, 5e1, 1e2, 1e3]
    for i, positions in enumerate(
        np.loadtxt(f).reshape(ITERATIONS // SAMPLING_RATE, N, 3)
        for f in sorted([f for f in Path().iterdir() if f.is_file() and "position" in f.name])
    ):
        mean_distances = np.empty(pos_size * pw_size)
        for j, x in enumerate(positions[CUTOFF:]):
            mean_distances[j * pw_size : (j + 1) * pw_size] = distance.pdist(x, "euclidean")
        sns.kdeplot(data=mean_distances, ax=ax, bw_adjust=0.2, linewidth=2, alpha=0.5, label=f"$T={temps[i]}$")
    ax.set(xlabel="Interparticle distance (a.u.)", ylabel="Pair correlation function", xlim=(0, 150))
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig("interparticle_distance.png")


def plot_number_density() -> None:
    """Plot the number density (y) against the distance (x)."""
    _, ax = plt.subplots(figsize=(12, 9))
    temps = [1e0, 1e1, 5e1, 1e2, 1e3]
    for (i, number_density), positions in zip(
        enumerate(
            np.loadtxt(f) / (ITERATIONS / SAMPLING_RATE - CUTOFF)
            for f in sorted([f for f in Path().iterdir() if f.is_file() and "density" in f.name and "out" in f.name])
        ),
        (np.loadtxt(f) for f in sorted([f for f in Path().iterdir() if f.is_file() and "position" in f.name])),
    ):
        distances = np.linalg.norm(positions[N * CUTOFF :, :2], axis=1)
        r_bin_length = 1.1 * (np.max(distances) - np.min(distances)) / BINS_X
        r = np.array([(r + 0.5) * r_bin_length for r in range(BINS_X)])
        n_norm = N * number_density / (2 * np.pi * r * np.trapz(number_density, dx=r_bin_length))
        ax.plot(r, n_norm, alpha=0.5, label=f"$T={temps[i]}$")
    ax.set(xlabel="Distance, $r$ (a.u.)", ylabel="Number density, $n(r)$ (a.u.)", xlim=(0, 4))
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig("number_density.png")


def plot_radial_distribution() -> None:
    """Plot the radial distribution (y) against the distance (x)."""
    _, ax = plt.subplots(figsize=(12, 9))
    temps = [1e0, 1e1, 5e1, 1e2, 1e3]
    for (i, distance), positions in zip(
        enumerate(
            np.loadtxt(f) / (ITERATIONS / SAMPLING_RATE - CUTOFF)
            for f in sorted([f for f in Path().iterdir() if f.is_file() and "density" in f.name and "out" in f.name])
        ),
        (np.loadtxt(f) for f in sorted([f for f in Path().iterdir() if f.is_file() and "position" in f.name])),
    ):
        distances = np.linalg.norm(positions[N * CUTOFF :, :2], axis=1)
        r_bin_length = 1.1 * (np.max(distances) - np.min(distances)) / BINS_X
        r = np.array([(r + 0.5) * r_bin_length for r in range(BINS_X)])
        n_norm = distance / np.trapz(distance, dx=r_bin_length)
        ax.plot(r, n_norm, alpha=0.5, label=f"$T={temps[i]}$")
    ax.set(xlabel="Distance", ylabel="Radial distribution")
    ax.legend(loc="upper right")
    plt.savefig("radial_distribution.png")


def plot_densities() -> None:
    """Plot numerous number densities from many runs of the Classical Metropolis algorithm."""
    plt.rcParams.update({"font.size": 20})
    _, ax = plt.subplots(figsize=(12, 9))
    temps = [1e-3, 2.5e-1, 5e-1, 1e0, 1e1]
    for (i, density), positions in zip(
        enumerate(
            np.loadtxt(f).reshape(BINS_X, BINS_Y) / (ITERATIONS / SAMPLING_RATE - CUTOFF)
            for f in sorted([f for f in Path().iterdir() if f.is_file() and "density" in f.name and "out" in f.name])
        ),
        (
            np.loadtxt(f).reshape(ITERATIONS // SAMPLING_RATE, N, 3)
            for f in sorted([f for f in Path().iterdir() if f.is_file() and "position" in f.name])
        ),
    ):
        x_bin_length = 1.1 * (np.max(positions[CUTOFF:, :, 0]) - np.min(positions[CUTOFF:, :, 0])) / BINS_X
        x = [x * x_bin_length * BINS_X for x in range(-BINS_X // 2 + 1, BINS_X // 2 + 1)]
        ax.plot(x, density[:, BINS_Y // 2], label=f"$T={temps[i]}$")
    ax.set(xlabel="$x$ position", ylabel="Number density (a.u.)", xlim=(-750, 750))
    ax.legend(loc="upper right")
    plt.savefig("densities.png")


if __name__ == "__main__":
    positions, energies = read_simulation_data()

    plot_energies_iterations(energies)
    # plot_many_energies_iterations()
    plot_density(positions)
    # plot_pair_density(positions)
    plot_snapshots(positions, -1)
