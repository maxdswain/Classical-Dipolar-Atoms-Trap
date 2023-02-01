#include <math.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_randist.h>
#include <toml.h>

#include "shared.h"

#define BAR_WIDTH 70
#define BOLTZMANN 1  // Boltzmann constant in defined systems of units based on values in the input file

typedef struct {
    int m, n, l;
    double *data;
} Double3D;

void read_config(int *N, int *ITERATIONS, double *T, double *M, double *SIGMA, double *DIPOLE_MOMENT,
                 double *DIPOLE_UNIT_VECTOR, double *FREQUENCY_Z, double *FREQUENCY_TRANSVERSE,
                 double *WALL_REPULSION_COEFFICIENT, int *WALL_REPULSION_ORDER, int *SAMPLING_RATE, int *READ_CONFIG,
                 int *SEED);
void position_random_generation(double **positions, int N, double T, double M, double FREQUENCY_Z,
                                double FREQUENCY_TRANSVERSE, int SEED);
double magnitude(double *a, int D);
double dot_product(double *a, double *b, int D);
void free_2D_array(double **array, int size);
double calculate_energy(double **positions, double *position, int N, double M, int index, double DIPOLE_MOMENT,
                        double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE,
                        double WALL_REPULSION_COEFFICIENT, int WALL_REPULSION_ORDER);
double calculate_total_energy(double **positions, int N, double M, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR,
                              double FREQUENCY_Z, double FREQUENCY_TRANSVERSE, double WALL_REPULSION_COEFFICIENT,
                              int WALL_REPULSION_ORDER);
double *metropolis_hastings(double **positions, int ITERATIONS, int N, double M, double T, double SIGMA,
                            double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z,
                            double FREQUENCY_TRANSVERSE, double WALL_REPULSION_COEFFICIENT, int WALL_REPULSION_ORDER,
                            int SAMPLING_RATE, int SEED);
void read_configuration(double **positions, int N);
void export_positions(Double3D *positions_saved);
void progress_bar(double progress, double time_taken);

int main(int argc, char **argv) {
    // Declare then read constants from an input file and print them out
    int N, ITERATIONS, SAMPLING_RATE, WALL_REPULSION_ORDER, READ_CONFIG, SEED;
    double T, M, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR[3], FREQUENCY_Z, FREQUENCY_TRANSVERSE,
        WALL_REPULSION_COEFFICIENT;

    read_config(&N, &ITERATIONS, &T, &M, &SIGMA, &DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, &FREQUENCY_Z,
                &FREQUENCY_TRANSVERSE, &WALL_REPULSION_COEFFICIENT, &WALL_REPULSION_ORDER, &SAMPLING_RATE, &READ_CONFIG,
                &SEED);
    printf(
        "Current variables set in input file:\nN: %d\niterations: %d\ntemperature: %f\nmass: %f\nsigma: "
        "%e\ndipole moment magnitude: %f\ndipole unit vector: %f %f %f\nfrequency_z: %f\nfrequency_transverse: "
        "%f\nwall repulsion coefficient: %f\nwall repulsion order: %d\n\n",
        N, ITERATIONS, T, M, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR[0], DIPOLE_UNIT_VECTOR[1], DIPOLE_UNIT_VECTOR[2],
        FREQUENCY_Z, FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT, WALL_REPULSION_ORDER);

    // Run simulation using values read from input file then calculate total energy of the last configuration
    double **positions = malloc(N * sizeof(**positions));
    for (int i = 0; i < N; i++) {
        positions[i] = malloc(3 * sizeof(*positions));
    }
    if (READ_CONFIG) {
        read_configuration(positions, N);
    } else {
        position_random_generation(positions, N, T, M, FREQUENCY_Z, FREQUENCY_TRANSVERSE, SEED);
    }
    double *energies_saved = metropolis_hastings(positions, ITERATIONS, N, M, T, SIGMA, DIPOLE_MOMENT,
                                                 DIPOLE_UNIT_VECTOR, FREQUENCY_Z, FREQUENCY_TRANSVERSE,
                                                 WALL_REPULSION_COEFFICIENT, WALL_REPULSION_ORDER, SAMPLING_RATE, SEED);

    printf("Energy of last sampled configuration: %e\n", energies_saved[ITERATIONS / SAMPLING_RATE - 1]);
    export_1D_array("energy_data.out", energies_saved, ITERATIONS / SAMPLING_RATE);
    free_2D_array(positions, N);
}

void read_config(int *N, int *ITERATIONS, double *T, double *M, double *SIGMA, double *DIPOLE_MOMENT,
                 double *DIPOLE_UNIT_VECTOR, double *FREQUENCY_Z, double *FREQUENCY_TRANSVERSE,
                 double *WALL_REPULSION_COEFFICIENT, int *WALL_REPULSION_ORDER, int *SAMPLING_RATE, int *READ_CONFIG,
                 int *SEED) {
    FILE *fp = fopen("input.toml", "r");  // Read and parse toml file
    if (!fp) {
        error("Cannot open input.toml", strerror(errno));
    }
    char errbuf[200];
    toml_table_t *conf = toml_parse_file(fp, errbuf, sizeof(errbuf));
    fclose(fp);
    if (!conf) {
        error("Cannot parse", errbuf);
    }
    toml_table_t *properties = toml_table_in(conf, "simulation_properties");  // Traverse to a table
    toml_datum_t particles = toml_int_in(properties, "particles");            // Extract values from table
    *N = particles.u.i;
    toml_datum_t repetitions = toml_int_in(properties, "repetitions");
    *ITERATIONS = repetitions.u.i;
    toml_datum_t temperature = toml_double_in(properties, "temperature");
    *T = temperature.u.d;
    toml_datum_t mass = toml_double_in(properties, "mass");
    *M = mass.u.d;
    toml_datum_t trial_sigma = toml_double_in(properties, "sigma");
    *SIGMA = trial_sigma.u.d;
    toml_datum_t dipole_moment_magnitude = toml_double_in(properties, "dipole_moment");
    *DIPOLE_MOMENT = dipole_moment_magnitude.u.d;
    toml_array_t *dipole_unit_vector_array = toml_array_in(properties, "dipole_unit_vector");
    for (int i = 0; i < 3; i++) {
        toml_datum_t toml_dipole_unit_vector = toml_double_at(dipole_unit_vector_array, i);
        DIPOLE_UNIT_VECTOR[i] = toml_dipole_unit_vector.u.d;
    }
    toml_datum_t trapping_frequency_z = toml_double_in(properties, "trapping_frequency_z");
    *FREQUENCY_Z = trapping_frequency_z.u.d;
    toml_datum_t trapping_frequency_transverse = toml_double_in(properties, "trapping_frequency_transverse");
    *FREQUENCY_TRANSVERSE = trapping_frequency_transverse.u.d;
    toml_datum_t wall_repulsion_coefficient = toml_double_in(properties, "wall_repulsion_coefficient");
    *WALL_REPULSION_COEFFICIENT = wall_repulsion_coefficient.u.d;
    toml_datum_t wall_repulsion_order = toml_int_in(properties, "order_repulsive_wall");
    *WALL_REPULSION_ORDER = wall_repulsion_order.u.i;
    toml_datum_t data_sampling_rate = toml_int_in(properties, "sampling_rate");
    *SAMPLING_RATE = data_sampling_rate.u.i;
    toml_datum_t read_configuration = toml_bool_in(properties, "monte_carlo_read_configuration");
    *READ_CONFIG = read_configuration.u.b;
    toml_datum_t seed = toml_int_in(properties, "seed");
    *SEED = seed.u.i;
    toml_free(conf);  // Free memory
}

void position_random_generation(double **positions, int N, double T, double M, double FREQUENCY_Z,
                                double FREQUENCY_TRANSVERSE, int SEED) {
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);  // Generator type
    gsl_rng_set(r, SEED);                         // Set seed

    // Randomly generate N x 3 array from uniform distribution
    const double r_xy = sqrt(6 * BOLTZMANN * T / (M * FREQUENCY_TRANSVERSE * FREQUENCY_TRANSVERSE));
    const double r_z = sqrt(6 * BOLTZMANN * T / (M * FREQUENCY_Z * FREQUENCY_Z));

    for (int i = 0; i < N; i++) {
        positions[i][0] = gsl_ran_flat(r, -r_xy, r_xy);
        positions[i][1] = gsl_ran_flat(r, -r_xy, r_xy);
        positions[i][2] = gsl_ran_flat(r, -r_z, r_z);
    }
    gsl_rng_free(r);
}

double magnitude(double *a, int D) {
    double result = 0;

    for (int i = 0; i < D; i++) {
        result += a[i] * a[i];
    }
    return sqrt(result);
}

double dot_product(double *a, double *b, int D) {
    double result = 0;

    for (int i = 0; i < D; i++) {
        result += a[i] * b[i];
    }
    return result;
}

void free_2D_array(double **array, int size) {
    for (int i = 0; i < size; i++) {
        free(array[i]);
    }
    free(array);
}

// Function that calculates the potential energy of a given single atom
double calculate_energy(double **positions, double *position, int N, double M, int index, double DIPOLE_MOMENT,
                        double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE,
                        double WALL_REPULSION_COEFFICIENT, int WALL_REPULSION_ORDER) {
    const double trapping_potential =
        0.5 * M *
        (FREQUENCY_TRANSVERSE * FREQUENCY_TRANSVERSE * (position[0] * position[0] + position[1] * position[1]) +
         FREQUENCY_Z * FREQUENCY_Z * position[2] * position[2]);
    double dipole_dipole_interaction = 0, wall_repulsion = 0;
    double displacement[3], distance, vector_term;

    for (int i = 0; i < N; i++) {
        if (i != index) {
            for (int j = 0; j < 3; j++) {
                displacement[j] = position[j] - positions[i][j];
            }
            distance = magnitude(displacement, 3);
            vector_term = dot_product(displacement, DIPOLE_UNIT_VECTOR, 3);
            dipole_dipole_interaction += (distance * distance - 3 * vector_term * vector_term) / pow(distance, 5);
            wall_repulsion += WALL_REPULSION_COEFFICIENT / pow(distance, WALL_REPULSION_ORDER);
        }
    }
    return trapping_potential + (DIPOLE_MOMENT * DIPOLE_MOMENT * dipole_dipole_interaction) + wall_repulsion;
}

// Function that calculates the total potential energy of a given configuration
double calculate_total_energy(double **positions, int N, double M, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR,
                              double FREQUENCY_Z, double FREQUENCY_TRANSVERSE, double WALL_REPULSION_COEFFICIENT,
                              int WALL_REPULSION_ORDER) {
    double total_dipole_dipole_interaction = 0, total_trapping_potential = 0, total_wall_repulsion = 0;
    double displacement[3], distance, vector_term;

    for (int i = 0; i < N; i++) {
        total_trapping_potential += 0.5 * M *
                                    (FREQUENCY_TRANSVERSE * FREQUENCY_TRANSVERSE *
                                         (positions[i][0] * positions[i][0] + positions[i][1] * positions[i][1]) +
                                     FREQUENCY_Z * FREQUENCY_Z * positions[i][2] * positions[i][2]);
        for (int j = 0; j < i; j++) {
            for (int k = 0; k < 3; k++) {
                displacement[k] = positions[i][k] - positions[j][k];
            }
            distance = magnitude(displacement, 3);
            vector_term = dot_product(displacement, DIPOLE_UNIT_VECTOR, 3);
            total_dipole_dipole_interaction += (distance * distance - 3 * vector_term * vector_term) / pow(distance, 5);
            total_wall_repulsion += WALL_REPULSION_COEFFICIENT / pow(distance, WALL_REPULSION_ORDER);
        }
    }
    return total_trapping_potential + (DIPOLE_MOMENT * DIPOLE_MOMENT * total_dipole_dipole_interaction) +
           total_wall_repulsion;
}

double *metropolis_hastings(double **positions, int ITERATIONS, int N, double M, double T, double SIGMA,
                            double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z,
                            double FREQUENCY_TRANSVERSE, double WALL_REPULSION_COEFFICIENT, int WALL_REPULSION_ORDER,
                            int SAMPLING_RATE, int SEED) {
    clock_t begin = clock();
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);  // Generator type
    gsl_rng_set(r, SEED);                         // Set seed

    int accepted = 0;
    double trial_positions[3], energy_difference, energy_previous;
    double *energies_saved = malloc(ITERATIONS / SAMPLING_RATE * sizeof(*energies_saved));
    Double3D positions_saved = {ITERATIONS / SAMPLING_RATE, N, 3};
    positions_saved.data =
        malloc(positions_saved.m * positions_saved.n * positions_saved.l * sizeof(*positions_saved.data));
    /* Main code of the Metropolis-Hasting algorithm. The previous energy is calculated then trial positions
    are generated for index particle, the difference in energy between the previous energy and the energy
    with the trial positions is calculated and if the energy different is less than or equal to 0 the trial
    move is accepted. If the energy difference is greater than or equal to a uniform randomly generated
    number between 0 and 1 then the trial move is also accepted, otherwise the move is rejected */
    for (int i = 0; i < ITERATIONS; i++) {
        for (int index = 0; index < N; index++) {
            energy_previous =
                calculate_energy(positions, positions[index], N, M, index, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR,
                                 FREQUENCY_Z, FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT, WALL_REPULSION_ORDER);
            for (int j = 0; j < 3; j++) {
                trial_positions[j] = positions[index][j] + gsl_ran_gaussian_ziggurat(r, SIGMA);
            }
            energy_difference =
                calculate_energy(positions, trial_positions, N, M, index, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR,
                                 FREQUENCY_Z, FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT, WALL_REPULSION_ORDER) -
                energy_previous;
            if (energy_difference <= 0) {
                accepted++;
                for (int j = 0; j < 3; j++) {
                    positions[index][j] = trial_positions[j];
                }
            } else if (gsl_rng_uniform(r) <= exp(-energy_difference / (BOLTZMANN * T))) {
                accepted++;
                for (int j = 0; j < 3; j++) {
                    positions[index][j] = trial_positions[j];
                }
            }
            /* Every SAMPLING_RATE iteration the positions of every atom are saved and
            the total energy of the configuration is calculated then saved */
            if (i % SAMPLING_RATE == 0) {
                energies_saved[i / SAMPLING_RATE] =
                    calculate_total_energy(positions, N, M, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, FREQUENCY_Z,
                                           FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT, WALL_REPULSION_ORDER);
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < 3; k++) {
                        positions_saved.data[i / SAMPLING_RATE * (positions_saved.n * positions_saved.l) +
                                             j * positions_saved.l + k] = positions[j][k];
                    }
                }
            }
        }
        double time_taken = (double)(clock() - begin) / CLOCKS_PER_SEC;
        progress_bar(((double)i + 1) / (double)ITERATIONS, time_taken);
    }
    gsl_rng_free(r);
    export_positions(&positions_saved);
    double percent_accepted = 100 * (double)accepted / (double)(ITERATIONS * N);
    printf("\n\npercent accepted: %.2f%%\nnumber accepted: %d\n\n\n", percent_accepted, accepted);
    return energies_saved;
}

void read_configuration(double **positions, int N) {
    FILE *fp = fopen("configuration.in", "r");
    if (!fp) {
        error("Cannot find configuration.in", strerror(errno));
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            fscanf(fp, "%lf", &positions[i][j]);
        }
    }
    fclose(fp);
}

void export_positions(Double3D *positions_saved) {
    FILE *fp = fopen("position_data.out", "w");
    for (int i = 0; i < positions_saved->m; i++) {
        for (int j = 0; j < positions_saved->n; j++) {
            fprintf(fp, "%f %f %f\n",
                    positions_saved->data[i * (positions_saved->n * positions_saved->l) + j * positions_saved->l],
                    positions_saved->data[i * (positions_saved->n * positions_saved->l) + j * positions_saved->l + 1],
                    positions_saved->data[i * (positions_saved->n * positions_saved->l) + j * positions_saved->l + 2]);
        }
    }
    fclose(fp);
    free(positions_saved->data);
}

void progress_bar(double progress, double time_taken) {
    int pos = BAR_WIDTH * progress;

    printf("\r|");
    for (int i = 0; i < BAR_WIDTH; i++) {
        if (i < pos) {
            printf("█");
        } else if (i == pos) {
            printf(">");
        } else {
            printf("-");
        }
    }
    if (progress == 1.0) {  // If complete, show time taken, otherwise don't
        printf("| %.2f %% Total time taken: %.2fs", progress * 100.0, time_taken);
    } else {
        printf("| %.2f %%", progress * 100.0);
    }
    fflush(stdout);
}
