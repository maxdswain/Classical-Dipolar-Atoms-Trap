#!/usr/bin/env python3

from collections.abc import Callable
from itertools import combinations
import os

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance
import seaborn as sns

# Read needed constants from input file
with open("input.toml", "rb") as f:
    input = tomllib.load(f)

N = input["simulation_properties"]["particles"]
ITERATIONS = input["simulation_properties"]["repetitions"]
T = input["simulation_properties"]["temperature"]
SAMPLING_RATE = input["simulation_properties"]["sampling_rate"]
ORDER_REPULSIVE_WALL = input["simulation_properties"]["order_repulsive_wall"]
WALL_REPULSION_COEFFICIENT = input["simulation_properties"]["wall_repulsion_coefficient"]
DIPOLE_MOMENT = np.linalg.norm(input["simulation_properties"]["dipole_vector"])
BINS_X = input["simulation_properties"]["bins_x"]
BINS_Y = input["simulation_properties"]["bins_y"]
BINS_Z = input["simulation_properties"]["bins_z"]
CUTOFF = input["simulation_properties"]["cutoff"]

plt.rcParams.update({"font.size": 20})
round_to_n = lambda x, n: x if x == 0 else round(x, -int(np.floor(np.log10(np.abs(x)))) + (n - 1))


def read_simulation_data() -> tuple[np.ndarray[np.float64], np.ndarray[np.float64]]:
    positions = np.loadtxt("position_data.out").reshape(ITERATIONS // SAMPLING_RATE, N, 3)
    energies = np.loadtxt("energy_data.out")
    return positions, energies


def boltzmann_distribution(energies: np.ndarray[np.float64]) -> Callable[[float], float]:
    BOLTZMANN_CONSTANT = 1  # Boltzmann constant in defined systems of units based on values in input.toml
    Z = sum([np.exp(-energy / (BOLTZMANN_CONSTANT * T)) for energy in energies])
    return lambda energy: len(energies) * np.exp(-energy / (BOLTZMANN_CONSTANT * T)) / Z


def plot_energy_histogram(energies: np.ndarray[np.float64]) -> None:
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


def plot_positions_iterations(positions: np.ndarray[np.float64], component: int) -> None:
    _, ax = plt.subplots()
    for i in range(N):
        ax.plot(np.linspace(SAMPLING_RATE, ITERATIONS, ITERATIONS // SAMPLING_RATE), positions[:, i, component])
    ax.set(xlabel="Iterations", ylabel=f"${['x', 'y', 'z'][component]}$ Positions")
    plt.show()


def plot_energies_iterations(energies: np.ndarray[np.float64]) -> None:
    _, ax = plt.subplots(figsize=(15, 9))
    ax.plot(np.linspace(0, ITERATIONS, len(energies)), energies)
    ax.set(xlabel="Iterations", ylabel="Energy", ylim=(0, 5 * round_to_n(energies[-1], 1)))
    plt.savefig("energies_iterations.png")


def plot_many_energies_iterations() -> None:
    _, ax = plt.subplots(figsize=(15, 9))
    energy_data = [np.loadtxt(f) for f in os.listdir(".") if os.path.isfile(f) and "energy" in f]
    for energies in energy_data:
        ax.plot(np.linspace(0, ITERATIONS, len(energies)), energies)
    ax.set(xlabel="Iterations", ylabel="Energy", xlim=(0, ITERATIONS))
    plt.ylim(1.117 * round_to_n(energy_data[0][-1], 1), 1.124 * round_to_n(energy_data[0][-1], 1))
    plt.savefig("many_energies_iterations.png")


def plot_snapshots(positions: np.ndarray[np.float64], iteration: int) -> None:
    plt.rcParams.update({"font.size": 14})
    _, axs = plt.subplots(1, 3, figsize=(15, 4))
    for i, (coord, x) in enumerate(zip(combinations("xyz", 2), combinations(range(3), 2))):
        axs[i].scatter(positions[iteration, :, x[0]], positions[iteration, :, x[1]], c="black")
        axs[i].set(xlabel=f"${coord[0]}$ Positions", ylabel=f"${coord[1]}$ Positions")
    plt.tight_layout()
    plt.savefig("snapshots.png")


def plot_error() -> None:
    _, ax = plt.subplots()
    errors = np.loadtxt("error_data.out")
    ax.plot(np.linspace(1, len(errors), len(errors), dtype="int64"), errors)
    ax.set(xlabel="BTN", ylabel="Standard Error")
    plt.show()


def plot_potential() -> None:
    _, ax = plt.subplots(figsize=(12, 9))
    r = np.linspace(0.9, 3, 101)  # Change ranges of values to ones relevant to system length scales
    repulsive_wall = WALL_REPULSION_COEFFICIENT * r**-ORDER_REPULSIVE_WALL
    dd_interaction = -2 * DIPOLE_MOMENT**2 * r**-3
    ax.plot(r, dd_interaction + repulsive_wall, label="Combined potential")
    ax.plot(r, repulsive_wall, label=rf"$\frac{{c_{6}}}{{r^{{{ORDER_REPULSIVE_WALL}}}}}$ potential")
    ax.plot(r, dd_interaction, label="Dipole-dipole interaction")
    ax.set(
        xlabel="Difference in position $r$",
        ylabel="Potential $U(r)$",
    )
    ax.legend(loc="upper right")
    plt.savefig("potential_components.png")


def plot_potentials() -> None:
    _, ax = plt.subplots()
    r = np.linspace(1, 6, 101)  # Change ranges of values to ones relevant to system length scales
    for x in [WALL_REPULSION_COEFFICIENT * 10**i for i in range(-2, 3)]:
        potential = x * r**-ORDER_REPULSIVE_WALL - 2 * DIPOLE_MOMENT**2 * r**-3
        ax.plot(r, potential, label=f"Potential for $c_{{6}}={x}$")
    ax.set(
        xlabel="Difference in Position $r$",
        ylabel="Potential $V(r)$",
        yscale="symlog",
    )
    ax.legend(loc="upper right")
    plt.show()


def plot_density(positions: np.ndarray[np.float64]) -> None:
    plt.rcParams.update({"font.size": 24})
    fig, ax = plt.subplots(figsize=(15, 12))
    density = np.sum(
        [
            np.loadtxt(f).reshape(BINS_X, BINS_Y)
            for f in os.listdir(".")
            if os.path.isfile(f) and "density" in f and "out" in f and "pair" not in f
        ],
        axis=0,
    )
    density /= ITERATIONS / SAMPLING_RATE - CUTOFF
    x_bin_length = 1.1 * (np.max(positions[CUTOFF:, :, 0]) - np.min(positions[CUTOFF:, :, 0])) / BINS_X
    y_bin_length = 1.1 * (np.max(positions[CUTOFF:, :, 1]) - np.min(positions[CUTOFF:, :, 1])) / BINS_Y
    cs = ax.contourf(density, 40, cmap="inferno")
    x_axes, y_axes = len(ax.get_xticklabels()), len(ax.get_yticklabels())
    # Note that positions labels are shifted to 0 to 2x from -x to x
    ax.set(
        xlabel="$x$ position",
        ylabel="$y$ position",
        xticklabels=[round_to_n(x * x_bin_length * BINS_X / x_axes, 4) for x in range(x_axes)],
        yticklabels=[round_to_n(y * y_bin_length * BINS_Y / y_axes, 4) for y in range(y_axes)],
    )
    fig.colorbar(cs).ax.set_ylabel("Number density $(a.u.)$", rotation=270, labelpad=30)
    plt.savefig("density_contour.png")


def plot_pair_density(positions: np.ndarray[np.float64]) -> None:
    plt.rcParams.update({"font.size": 24})
    fig, ax = plt.subplots(figsize=(15, 12))
    pair_density = np.sum(
        [
            np.loadtxt(f).reshape(BINS_Z, BINS_Z)
            for f in os.listdir(".")
            if os.path.isfile(f) and "pair" in f and "out" in f
        ],
        axis=0,
    )
    pair_density /= ITERATIONS / SAMPLING_RATE - CUTOFF
    z_bin_length = 1.1 * (np.max(positions[CUTOFF:, :, 2]) - np.min(positions[CUTOFF:, :, 2])) / BINS_Z
    cs = ax.contourf(pair_density, 40, cmap="inferno")
    x_axes, y_axes = len(ax.get_xticklabels()), len(ax.get_yticklabels())
    # Note that positions labels are shifted to 0 to 2x from -x to x
    ax.set(
        xlabel="$z_{1}$ position",
        ylabel="$z_{2}$ position",
        xticklabels=[round_to_n(z * z_bin_length * BINS_Z / x_axes, 4) for z in range(x_axes)],
        yticklabels=[round_to_n(z * z_bin_length * BINS_Z / y_axes, 4) for z in range(y_axes)],
    )
    fig.colorbar(cs).ax.set_ylabel("Mean number of pairs", rotation=270, labelpad=30)
    plt.savefig("pair_density_contour.png")


def plot_interparticle_distance() -> None:
    _, ax = plt.subplots(figsize=(12, 9))
    pw_size, pos_size = int(0.5 * N * (N - 1)), int(ITERATIONS / SAMPLING_RATE - CUTOFF)
    temps = [1e-3, 2.5e-1, 5e-1, 1e0, 1e1]
    for i, positions in enumerate(
        map(
            lambda f: np.loadtxt(f).reshape(ITERATIONS // SAMPLING_RATE, N, 3),
            sorted([f for f in os.listdir(".") if os.path.isfile(f) and "position" in f]),
        )
    ):
        mean_distances = np.empty(pos_size * pw_size)
        for j, x in enumerate(positions[CUTOFF:]):
            mean_distances[j * pw_size : (j + 1) * pw_size] = distance.pdist(x, "euclidean")
        sns.kdeplot(data=mean_distances, ax=ax, bw_adjust=0.2, linewidth=2, alpha=0.5, label=f"$T={temps[i]}$")
    ax.set(xlabel="Interparticle distance", ylabel="Density", xlim=(0, 20))
    ax.legend(loc="upper right")
    plt.savefig("interparticle_distance.png")


if __name__ == "__main__":
    positions, energies = read_simulation_data()

    plot_energies_iterations(energies)
    # plot_many_energies_iterations()
    # plot_density(positions)
    # plot_pair_density(positions)
    plot_snapshots(positions, -1)
