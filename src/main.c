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

void read_config(int *N, int *ITERATIONS, double *T, double *M, double *SIGMA, double *DIPOLE_VECTOR,
                 double *FREQUENCY_Z, double *FREQUENCY_TRANSVERSE, double *WALL_COEFFICIENT, int *WALL_ORDER,
                 int *SAMPLING_RATE, int *READ_CONFIG, int *SEED);
void position_random_generation(double **positions, int N, double T, double M, double FREQUENCY_Z,
                                double FREQUENCY_TRANSVERSE, int SEED);
double magnitude(double *a, int D);
double dot_product(double *a, double *b, int D);
double calculate_energy(double **positions, double *position, int N, double M, int index, double *DIPOLE_VECTOR,
                        double FREQUENCY_Z, double FREQUENCY_TRANSVERSE, double WALL_COEFFICIENT, int WALL_ORDER);
double calculate_total_energy(double **positions, int N, double M, double *DIPOLE_VECTOR, double FREQUENCY_Z,
                              double FREQUENCY_TRANSVERSE, double WALL_COEFFICIENT, int WALL_ORDER);
double *metropolis_hastings(double **positions, int ITERATIONS, int N, double M, double T, double SIGMA,
                            double *DIPOLE_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE,
                            double WALL_COEFFICIENT, int WALL_ORDER, int SAMPLING_RATE, int SEED);
void export_positions(Double3D *pos_saved);
void progress_bar(double progress, double time_taken);

int main(int argc, char **argv) {
    // Declare then read constants from an input file and print them out
    int N, ITERATIONS, SAMPLING_RATE, WALL_ORDER, READ_CONFIG, SEED;
    double T, M, SIGMA, DIPOLE_VECTOR[3], FREQUENCY_Z, FREQUENCY_TRANSVERSE, WALL_COEFFICIENT;

    read_config(&N, &ITERATIONS, &T, &M, &SIGMA, DIPOLE_VECTOR, &FREQUENCY_Z, &FREQUENCY_TRANSVERSE, &WALL_COEFFICIENT,
                &WALL_ORDER, &SAMPLING_RATE, &READ_CONFIG, &SEED);
    printf(
        "Current variables set in input file:\nN: %d\niterations: %d\ntemperature: %f\nmass: %f\nsigma: "
        "%e\ndipole vector: %f %f %f\nfrequency_z: %f\nfrequency_transverse: "
        "%f\nwall repulsion coefficient: %f\nwall repulsion order: %d\n\n",
        N, ITERATIONS, T, M, SIGMA, DIPOLE_VECTOR[0], DIPOLE_VECTOR[1], DIPOLE_VECTOR[2], FREQUENCY_Z,
        FREQUENCY_TRANSVERSE, WALL_COEFFICIENT, WALL_ORDER);

    // Run simulation using values read from input file then calculate total energy of the last configuration
    double **positions = malloc(N * sizeof(**positions));
    for (int i = 0; i < N; i++) {
        positions[i] = malloc(3 * sizeof(*positions));
    }
    if (READ_CONFIG) {
        read_2D_array("configuration.in", positions, N);
    } else {
        position_random_generation(positions, N, T, M, FREQUENCY_Z, FREQUENCY_TRANSVERSE, SEED);
    }
    double *energies_saved =
        metropolis_hastings(positions, ITERATIONS, N, M, T, SIGMA, DIPOLE_VECTOR, FREQUENCY_Z, FREQUENCY_TRANSVERSE,
                            WALL_COEFFICIENT, WALL_ORDER, SAMPLING_RATE, SEED);

    printf("Energy of last sampled configuration: %e\n", energies_saved[ITERATIONS / SAMPLING_RATE - 1]);
    export_1D_array("energy_data", energies_saved, ITERATIONS / SAMPLING_RATE);
    free_2D_array(positions, N);
}

void read_config(int *N, int *ITERATIONS, double *T, double *M, double *SIGMA, double *DIPOLE_VECTOR,
                 double *FREQUENCY_Z, double *FREQUENCY_TRANSVERSE, double *WALL_COEFFICIENT, int *WALL_ORDER,
                 int *SAMPLING_RATE, int *READ_CONFIG, int *SEED) {
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
    toml_table_t *system = toml_table_in(conf, "system");          // Traverse to a table
    toml_table_t *metropolis = toml_table_in(conf, "metropolis");  // Traverse to a table
    toml_datum_t particles = toml_int_in(system, "particles");     // Extract values from table
    *N = particles.u.i;
    toml_datum_t repetitions = toml_int_in(metropolis, "iterations");
    *ITERATIONS = repetitions.u.i;
    toml_datum_t temperature = toml_double_in(system, "temperature");
    *T = temperature.u.d;
    toml_datum_t mass = toml_double_in(system, "mass");
    *M = mass.u.d;
    toml_datum_t trial_sigma = toml_double_in(metropolis, "sigma");
    *SIGMA = trial_sigma.u.d;
    toml_array_t *dipole_vector_array = toml_array_in(system, "dipole_vector");
    for (int i = 0; i < 3; i++) {
        toml_datum_t toml_dipole_vector = toml_double_at(dipole_vector_array, i);
        DIPOLE_VECTOR[i] = toml_dipole_vector.u.d;
    }
    toml_datum_t trapping_frequency_z = toml_double_in(system, "trapping_frequency_z");
    *FREQUENCY_Z = trapping_frequency_z.u.d;
    toml_datum_t trapping_frequency_transverse = toml_double_in(system, "trapping_frequency_transverse");
    *FREQUENCY_TRANSVERSE = trapping_frequency_transverse.u.d;
    toml_datum_t wall_repulsion_coefficient = toml_double_in(system, "wall_coefficient");
    *WALL_COEFFICIENT = wall_repulsion_coefficient.u.d;
    toml_datum_t wall_repulsion_order = toml_int_in(system, "wall_order");
    *WALL_ORDER = wall_repulsion_order.u.i;
    toml_datum_t data_sampling_rate = toml_int_in(metropolis, "sampling_rate");
    *SAMPLING_RATE = data_sampling_rate.u.i;
    toml_datum_t read_configuration = toml_bool_in(metropolis, "read_initial_configuration");
    *READ_CONFIG = read_configuration.u.b;
    toml_datum_t seed = toml_int_in(metropolis, "seed");
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

// Function that calculates the potential energy of a given single atom
double calculate_energy(double **positions, double *position, int N, double M, int index, double *DIPOLE_VECTOR,
                        double FREQUENCY_Z, double FREQUENCY_TRANSVERSE, double WALL_COEFFICIENT, int WALL_ORDER) {
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
            vector_term = dot_product(displacement, DIPOLE_VECTOR, 3);
            dipole_dipole_interaction += (distance * distance - 3 * vector_term * vector_term) / pow(distance, 5);
            wall_repulsion += 1 / pow(distance, WALL_ORDER);
        }
    }
    return trapping_potential + dipole_dipole_interaction + WALL_COEFFICIENT * wall_repulsion;
}

// Function that calculates the total potential energy of a given configuration
double calculate_total_energy(double **positions, int N, double M, double *DIPOLE_VECTOR, double FREQUENCY_Z,
                              double FREQUENCY_TRANSVERSE, double WALL_COEFFICIENT, int WALL_ORDER) {
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
            vector_term = dot_product(displacement, DIPOLE_VECTOR, 3);
            total_dipole_dipole_interaction += (distance * distance - 3 * vector_term * vector_term) / pow(distance, 5);
            total_wall_repulsion += 1 / pow(distance, WALL_ORDER);
        }
    }
    return total_trapping_potential + total_dipole_dipole_interaction + WALL_COEFFICIENT * total_wall_repulsion;
}

double *metropolis_hastings(double **positions, int ITERATIONS, int N, double M, double T, double SIGMA,
                            double *DIPOLE_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE,
                            double WALL_COEFFICIENT, int WALL_ORDER, int SAMPLING_RATE, int SEED) {
    clock_t begin = clock();
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);  // Generator type
    gsl_rng_set(r, SEED);                         // Set seed

    int accepted = 0;
    double trial_positions[3], energy_difference, energy_previous;
    double *energies_saved = malloc(ITERATIONS / SAMPLING_RATE * sizeof(*energies_saved));
    Double3D pos_saved = {ITERATIONS / SAMPLING_RATE, N, 3};
    pos_saved.data = malloc(pos_saved.m * pos_saved.n * pos_saved.l * sizeof(*pos_saved.data));
    /* Main code of the Metropolis-Hasting algorithm. The previous energy is calculated then trial positions
    are generated for index particle, the difference in energy between the previous energy and the energy
    with the trial positions is calculated and if the energy different is less than or equal to 0 the trial
    move is accepted. If the energy difference is greater than or equal to a uniform randomly generated
    number between 0 and 1 then the trial move is also accepted, otherwise the move is rejected */
    for (int i = 0; i < ITERATIONS; i++) {
        for (int index = 0; index < N; index++) {
            energy_previous = calculate_energy(positions, positions[index], N, M, index, DIPOLE_VECTOR, FREQUENCY_Z,
                                               FREQUENCY_TRANSVERSE, WALL_COEFFICIENT, WALL_ORDER);
            for (int j = 0; j < 3; j++) {
                trial_positions[j] = positions[index][j] + gsl_ran_gaussian_ziggurat(r, SIGMA);
            }
            energy_difference = calculate_energy(positions, trial_positions, N, M, index, DIPOLE_VECTOR, FREQUENCY_Z,
                                                 FREQUENCY_TRANSVERSE, WALL_COEFFICIENT, WALL_ORDER) -
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
                energies_saved[i / SAMPLING_RATE] = calculate_total_energy(
                    positions, N, M, DIPOLE_VECTOR, FREQUENCY_Z, FREQUENCY_TRANSVERSE, WALL_COEFFICIENT, WALL_ORDER);
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < 3; k++) {
                        pos_saved.data[i / SAMPLING_RATE * (pos_saved.n * pos_saved.l) + j * pos_saved.l + k] =
                            positions[j][k];
                    }
                }
            }
        }
        double time_taken = (double)(clock() - begin) / CLOCKS_PER_SEC;
        progress_bar(((double)i + 1) / (double)ITERATIONS, time_taken);
    }
    gsl_rng_free(r);
    export_positions(&pos_saved);
    double percent_accepted = 100 * (double)accepted / (double)(ITERATIONS * N);
    printf("\n\npercent accepted: %.2f%%\nnumber accepted: %d\n\n\n", percent_accepted, accepted);
    return energies_saved;
}

void export_positions(Double3D *pos_saved) {
    char *file_name = produce_file_name("position_data");
    FILE *fp = fopen(file_name, "w");
    for (int i = 0; i < pos_saved->m; i++) {
        for (int j = 0; j < pos_saved->n; j++) {
            fprintf(fp, "%f %f %f\n", pos_saved->data[i * (pos_saved->n * pos_saved->l) + j * pos_saved->l],
                    pos_saved->data[i * (pos_saved->n * pos_saved->l) + j * pos_saved->l + 1],
                    pos_saved->data[i * (pos_saved->n * pos_saved->l) + j * pos_saved->l + 2]);
        }
    }
    fclose(fp);
    free(file_name);
    free(pos_saved->data);
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
