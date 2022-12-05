#!/usr/bin/env python3

from collections.abc import Callable

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import matplotlib.pyplot as plt
import numpy as np

# Read needed constants from config file
with open("config.toml", "rb") as f:
    config = tomllib.load(f)

N = config["simulation_properties"]["particles"]
ITERATIONS = config["simulation_properties"]["repetitions"]
T = config["simulation_properties"]["temperature"]
SAMPLING_RATE = config["simulation_properties"]["sampling_rate"]
BTN = config["simulation_properties"]["blocking_transformation_number"]
CUTOFF = config["simulation_properties"]["cutoff"]


def read_simulation_data() -> tuple[np.ndarray[np.float64], np.ndarray[np.float64], np.ndarray[np.float64]]:
    positions = np.loadtxt("simulation_position_data.txt").reshape(ITERATIONS // SAMPLING_RATE, N, 3)
    energies = np.loadtxt("simulation_energy_data.txt")
    errors = np.loadtxt("simulation_error_data.txt")
    return positions, energies, errors


def boltzmann_distribution(energies: np.ndarray[np.float64]) -> Callable[[float], float]:
    BOLTZMANN_CONSTANT = 2617360049  # Boltzmann constant in defined systems of units based on values in config
    Z = sum([np.exp(-energy / (BOLTZMANN_CONSTANT * T)) for energy in energies])
    return lambda energy: len(energies) * np.exp(-energy / (BOLTZMANN_CONSTANT * T)) / Z


def plot_energy_histogram(energies: np.ndarray[np.float64]) -> None:
    _, ax = plt.subplots()
    sorted_energies = sorted(energies)
    distribution = boltzmann_distribution(sorted_energies)
    ax.hist(energies, bins=30, alpha=0.5)
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


def plot_energies_iterations(energies: np.ndarray[np.float64], error: float) -> None:
    _, ax = plt.subplots()
    ax.errorbar(
        np.linspace(CUTOFF * SAMPLING_RATE, ITERATIONS, len(energies)), energies, yerr=error, ecolor="r", fmt="k"
    )
    ax.set(xlabel="Iterations", ylabel="Mean energy")
    plt.show()


def plot_snapshot(positions: np.ndarray[np.float64], iteration: int) -> None:
    _, ax = plt.subplots()
    ax.scatter(positions[iteration, :, 0], positions[iteration, :, 1], c="black")
    ax.set(xlabel="$x$ Positions", ylabel="$y$ Positions")
    plt.show()


def plot_error(errors: np.ndarray[np.float64]) -> None:
    _, ax = plt.subplots()
    ax.plot(np.linspace(1, len(errors), len(errors), dtype="int64"), errors)
    ax.set(xlabel="BTN", ylabel="Standard Error")
    plt.show()


if __name__ == "__main__":
    positions, energies, errors = read_simulation_data()
    distances = np.linalg.norm(positions[-1], axis=1)

    plot_error(errors)
    plot_energy_histogram(energies)
    plot_energies_iterations(energies, errors[BTN - 1])
    plot_positions_iterations(positions, 0)
