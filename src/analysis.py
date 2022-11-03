#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import tomli

with open("config.toml", "rb") as f:
    config = tomli.load(f)

N = config["simulation_properties"]["particles"]
ITERATIONS = config["simulation_properties"]["repetitions"]
T = config["simulation_properties"]["temperature"]
M = config["simulation_properties"]["mass"]
DIPOLE_MOMENT = config["simulation_properties"]["dipole_moment"]
DIPOLE_UNIT_VECTOR = np.array(config["simulation_properties"]["dipole_unit_vector"])
FREQUENCY_Z = config["simulation_properties"]["trapping_frequency_z"]
FREQUENCY_TRANSVERSE = config["simulation_properties"]["trapping_frequency_transverse"]

def read_simulation_data() -> npt.NDArray[np.float64]:
    positions = np.empty(shape=(ITERATIONS, N, 3))
    with open("simulation_data.txt") as f:
        for i, line in enumerate(f):
            positions[int(i / N)][i % N] = [float(x) for x in line.split()]
    return positions

def calculate_energies(positions: npt.NDArray[np.float64]) -> list:
    energies = []
    for i in range(N):
        trapping_potential = 0.5 * M * FREQUENCY_Z**2 * ((FREQUENCY_TRANSVERSE / FREQUENCY_Z)**2 * np.sum(positions[-1][i][:2]**2) + positions[-1][i][2]**2)
        temp = 0
        for j in range(N):
            if j < i:
                displacement = positions[-1][i] - positions[-1][j]
                distance = np.linalg.norm(displacement)
                temp += distance**-6 + DIPOLE_MOMENT**2 * (distance**2 - 3 * np.dot(displacement, DIPOLE_UNIT_VECTOR)**2) / distance**5
        energies.append(trapping_potential + temp)
    return energies

if __name__ == "__main__":
    positions = read_simulation_data()
    energies = calculate_energies(positions)
    distances = np.linalg.norm(positions[-1], axis=1)

    _, ax = plt.subplots(1, 1)
    ax.hist(energies, bins=30, density=True, alpha=0.5)
    # ax.plot(PLACEHOLDER1, PLACEHOLDER2, "r-", lw=2, label="Boltzmann distribution")
    ax.set(xlabel="Energies", ylabel="Frequency")
    # ax.legend(loc="upper right")
    plt.show()

    print(positions[-1], sum(energies), distances)
