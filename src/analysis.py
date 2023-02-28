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

plt.rcParams.update({"font.size": 22})
round_to_n = lambda x, n: round(x, -int(np.floor(np.log10(np.abs(x)))) + (n - 1))


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
    ax.set(xlabel="Iterations", ylabel="Energy")
    plt.ylim(0, 5 * round_to_n(energies[-1], 1))
    plt.savefig("energies_iterations.png")


def plot_many_energies_iterations() -> None:
    _, ax = plt.subplots(figsize=(15, 9))
    energy_data = [np.loadtxt(f) for f in os.listdir(".") if os.path.isfile(f) and "energy" in f]
    for energies in energy_data:
        ax.plot(np.linspace(0, ITERATIONS, len(energies)), energies)
    ax.set(xlabel="Iterations", ylabel="Energy")
    plt.ylim(0.9 * round_to_n(energy_data[0][-1], 1), 1.5 * round_to_n(energy_data[0][-1], 1))
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
    _, ax = plt.subplots()
    r = np.linspace(0.1, 0.3, 101)  # Change ranges of values to ones relevant to system length scales
    repulsive_wall = WALL_REPULSION_COEFFICIENT * r**-ORDER_REPULSIVE_WALL
    dd_interaction = -2 * DIPOLE_MOMENT**2 * r**-3
    ax.plot(r, dd_interaction + repulsive_wall, label="Combined Potential")
    ax.plot(r, repulsive_wall, label=rf"$\frac{{1}}{{r^{{{ORDER_REPULSIVE_WALL}}}}}$ Potential")
    ax.plot(r, dd_interaction, label="Dipole-Diple Interaction")
    ax.set(
        xlabel="Difference in Position $r$",
        ylabel="Potential $V(r)$",
        title="Plot of Potentials Against the Difference in Distance Between 2 Atoms",
    )
    ax.legend(loc="upper right")
    plt.show()


def plot_potentials() -> None:
    _, ax = plt.subplots()
    r = np.linspace(1, 6, 101)  # Change ranges of values to ones relevant to system length scales
    for x in [WALL_REPULSION_COEFFICIENT * 10**i for i in range(-2, 3)]:
        potential = x * r**-ORDER_REPULSIVE_WALL - 2 * DIPOLE_MOMENT**2 * r**-3
        ax.plot(r, potential, label=f"Potential for $c_{{6}}={x}$")
    ax.set(
        xlabel="Difference in Position $r$",
        ylabel="Potential $V(r)$",
        title="Plot of Potentials Against the Difference in Distance Between 2 Atoms",
        yscale="symlog",
    )
    ax.legend(loc="upper right")
    plt.show()


def plot_density() -> None:
    fig, ax = plt.subplots(figsize=(15, 12))
    density = np.sum(
        [
            np.loadtxt(f).reshape(BINS_X, BINS_Y)
            for f in os.listdir(".")
            if os.path.isfile(f) and "density" in f and "out" in f and "pair" not in f
        ],
        axis=0,
    )
    cs = ax.contourf(density, 40, cmap="RdGy")
    ax.set(yticklabels=[], xticklabels=[])
    fig.colorbar(cs)
    plt.savefig("density_contour.png")


def plot_pair_density() -> None:
    fig, ax = plt.subplots(figsize=(15, 12))
    pair_density = np.sum(
        [
            np.loadtxt(f).reshape(BINS_Z, BINS_Z)
            for f in os.listdir(".")
            if os.path.isfile(f) and "pair" in f and "out" in f
        ],
        axis=0,
    )
    cs = ax.contourf(pair_density, 40, cmap="RdGy")
    ax.set(yticklabels=[], xticklabels=[])
    fig.colorbar(cs)
    plt.savefig("pair_density_contour.png")


if __name__ == "__main__":
    positions, energies = read_simulation_data()
    distances = np.linalg.norm(positions[-1], axis=1)

    plot_energies_iterations(energies)
    # plot_many_energies_iterations()
    # plot_density()
    # plot_pair_density()
    plot_snapshots(positions, -1)
