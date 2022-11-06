#!/usr/bin/env python3

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

with open("config.toml", "rb") as f:
    config = tomllib.load(f)

N = config["simulation_properties"]["particles"]
ITERATIONS = config["simulation_properties"]["repetitions"]

def read_simulation_data() -> tuple[npt.NDArray[np.float64], list[float]]:
    positions = np.empty(shape=(ITERATIONS, N, 3))
    with open("simulation_position_data.txt") as f:
        for i, line in enumerate(f):
            positions[int(i / N)][i % N] = [float(x) for x in line.split()]
    with open("simulation_energy_data.txt") as f:
        energies = [float(line) for line in f]
    return positions, energies

def plot_energy_histogram(energies: list[float]) -> None:
    _, ax = plt.subplots(1, 1)
    ax.hist(energies, bins=30, density=True, alpha=0.5)
    # ax.plot(PLACEHOLDER1, PLACEHOLDER2, "r-", lw=2, label="Boltzmann distribution")
    ax.set(xlabel="Energies", ylabel="Frequency")
    # ax.legend(loc="upper right")
    plt.show()

def plot_positions_iterations(positions: npt.NDArray[np.float64], component: int) -> None:
    _, ax = plt.subplots(1, 1)
    ax.scatter([i for i in range(1, ITERATIONS + 1) for j in range(N)], positions[..., component].reshape(ITERATIONS * N,), alpha=0.5)
    plt.show()

if __name__ == "__main__":
    positions, energies = read_simulation_data()
    distances = np.linalg.norm(positions[-1], axis=1)

    # plot_energy_histogram(energies)
    plot_positions_iterations(positions, 0)

    print(positions[-1][19][0], sum(energies))
