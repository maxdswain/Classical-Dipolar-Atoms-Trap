#!/usr/bin/env python3

from collections.abc import Callable

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import matplotlib.pyplot as plt
import numpy as np

with open("config.toml", "rb") as f:
    config = tomllib.load(f)

N = config["simulation_properties"]["particles"]
ITERATIONS = config["simulation_properties"]["repetitions"]
T = config["simulation_properties"]["temperature"]
SAMPLING_RATE = config["simulation_properties"]["sampling_rate"]

def read_simulation_data() -> tuple[np.ndarray[np.float64], list[float]]:
    positions = np.empty(shape=(ITERATIONS // SAMPLING_RATE, N, 3))
    with open("simulation_position_data.txt") as f:
        for i, line in enumerate(f):
            positions[int(i / N)][i % N] = [float(x) for x in line.split()]
    with open("simulation_energy_data.txt") as f:
        energies = [float(line) for line in f]
    return positions, energies[500:]

def boltzmann_distribution(energies: list[float]) -> Callable[[float], float]:
    BOLTZMANN_CONSTANT = 3.167e-6  # Boltzmann constant in Hartree units
    Z = sum([np.exp(-energy / (BOLTZMANN_CONSTANT * T)) for energy in energies])
    return lambda energy: len(energies) * np.exp(-energy / (BOLTZMANN_CONSTANT * T)) / Z

def plot_energy_histogram(energies: list[float]) -> None:
    _, ax = plt.subplots()
    sorted_energies = sorted(energies)
    distribution = boltzmann_distribution(sorted_energies)
    ax.hist(energies, bins=30, alpha=0.5)
    ax.plot(sorted_energies, [distribution(energy) for energy in sorted_energies], "r-", lw=2, label="Boltzmann distribution")
    ax.set(xlabel="Energy", ylabel="Frequency")
    ax.legend(loc="upper right")
    plt.show()

def plot_positions_iterations(positions: np.ndarray[np.float64], component: int) -> None:
    _, ax = plt.subplots()
    for i in range(N):
        ax.plot([SAMPLING_RATE * i for i in range(1, ITERATIONS // SAMPLING_RATE + 1)], positions[:, i, component])
    ax.set(xlabel="Iterations", ylabel=f"${['x', 'y', 'z'][component]}$ Positions")
    # plt.scatter([i for i in range(1, ITERATIONS + 1) for j in range(N)], positions[..., component].reshape(ITERATIONS * N,), alpha=0.5)
    plt.show()

if __name__ == "__main__":
    positions, energies = read_simulation_data()
    distances = np.linalg.norm(positions[-1], axis=1)

    plot_energy_histogram(energies)
    plot_positions_iterations(positions, 0)
