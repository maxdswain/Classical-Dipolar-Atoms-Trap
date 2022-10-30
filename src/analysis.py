#!/usr/bin/env python3

from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import tomli
from scipy.stats import maxwell

with open("config.toml", "rb") as f:
    config = tomli.load(f)

N = config["simulation_properties"]["particles"]
T = config["simulation_properties"]["temperature"]
m = config["simulation_properties"]["mass"]
dipole_moment = config["simulation_properties"]["magnetic_dipole_moment"]
frequency_z = config["simulation_properties"]["trapping_frequency_z"]
frequency_transverse = config["simulation_properties"]["trapping_frequency_transverse"]

def read_simulation_data():
    positions, counter, temp_array = [], 0, np.empty(shape=(N, 3))
    with open("simulation_data.txt") as f:
        for line in f:
            temp_array[counter] = [float(x) for x in line.split()]
            counter += 1
            if counter == N:
                positions.append(deepcopy(temp_array))
                counter, temp_array = 0, np.empty(shape=(N, 3))
    return positions

def calculate_energies(positions):
    mu_zero = 1.25663706212e-6
    temp = 0
    energies = []
    dipole_unit_vector = np.array([8 / np.sqrt(74), 3 / np.sqrt(74), 1 / np.sqrt(74)])
    for i in range(N):
        trapping_potential = 0.5 * m * frequency_z**2 * ((frequency_transverse / frequency_z)**2 * (positions[-1][i][0]**2 + positions[-1][i][1]**2) + positions[-1][i][2]**2)
        for j in range(N):
            if j < i:
                displacement = positions[-1][i] - positions[-1][j]
                distance = np.linalg.norm(displacement)
                vector_term = np.dot(displacement, dipole_unit_vector)
                temp += 1 / distance**12 + (mu_zero * dipole_moment * dipole_moment) * (distance**2 - 3 * vector_term**2) / (4 * np.pi * distance**5)
        energies.append(trapping_potential + temp)
        temp = 0
    return energies

if __name__ == "__main__":
    positions = read_simulation_data()
    energies = calculate_energies(positions)
    distances = np.linalg.norm(positions[-1], axis=1)

    _, ax = plt.subplots(1, 1)
    ax.hist(distances)
    # ax.plot()
    # ax.set(xlabel="Energy", ylabel="Frequency")
    # ax.legend(loc="upper right")
    plt.show()

    print(positions[0], sum(energies), distances)
