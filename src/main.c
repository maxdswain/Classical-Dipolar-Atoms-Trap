#include <math.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_randist.h>
#include <toml.h>

#define BAR_WIDTH 70
#define BOLTZMANN 2617360049  // Boltzmann constant in defined systems of units based on values in config

typedef struct {
    int m, n, l;
    double *data;
} Double3D;

void error(const char *msg, const char *errmsg);
void read_config(int *N, int *ITERATIONS, double *T, double *M, double *SIGMA, double *DIPOLE_MOMENT,
                 double *DIPOLE_UNIT_VECTOR, double *FREQUENCY_Z, double *FREQUENCY_TRANSVERSE,
                 double *WALL_REPULSION_COEFFICIENT, int *SAMPLING_RATE, int *BTN, int *CUTOFF);
double **position_random_generation(int N, double T, double M, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE);
double sum(double *a, int D);
double magnitude(double *a, int D);
double dot_product(double *a, double *b, int D);
double *slice(double *a, int from, int until);
void free_2D_array(double **array, int size);
double calculate_energy(double **positions, double *position, int N, double M, int index, double DIPOLE_MOMENT,
                        double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE,
                        double WALL_REPULSION_COEFFICIENT);
double calculate_total_energy(double **positions, int N, double M, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR,
                              double FREQUENCY_Z, double FREQUENCY_TRANSVERSE, double WALL_REPULSION_COEFFICIENT);
double *metropolis_hastings(double **positions, int ITERATIONS, int N, double M, double T, double SIGMA,
                            double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z,
                            double FREQUENCY_TRANSVERSE, double WALL_REPULSION_COEFFICIENT, int SAMPLING_RATE);
double *reblocking(double *energies_saved, int size, int BTN);
double calculate_error(double *mean_energies, int size, int N);
void export_positions(Double3D *positions_saved);
void export_1D_array(char *file_name, double *array, int size);
void progress_bar(double progress, double time_taken);

int main(int argc, char **argv) {
    // Declare then read constants from a config file and print them out
    int N, ITERATIONS, SAMPLING_RATE, BTN, CUTOFF;
    double T, M, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR[3], FREQUENCY_Z, FREQUENCY_TRANSVERSE,
        WALL_REPULSION_COEFFICIENT;

    read_config(&N, &ITERATIONS, &T, &M, &SIGMA, &DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, &FREQUENCY_Z,
                &FREQUENCY_TRANSVERSE, &WALL_REPULSION_COEFFICIENT, &SAMPLING_RATE, &BTN, &CUTOFF);
    printf(
        "Current variables set in config:\nN: %d\niterations: %d\ntemperature: %f\nmass: %f\nsigma: "
        "%e\ndipole moment magnitude: %f\ndipole unit vector: %f %f %f\nfrequency_z: %f\nfrequency_transverse: "
        "%f\nwall repulsion coefficient: %f\ncutoff: %d\n\n",
        N, ITERATIONS, T, M, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR[0], DIPOLE_UNIT_VECTOR[1], DIPOLE_UNIT_VECTOR[2],
        FREQUENCY_Z, FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT, CUTOFF);

    int MAX_BTN = 2;  // ((ITERATIONS / SAMPLING_RATE) - CUTOFF) / 2^MAX_BTN must be an integer
    while (((ITERATIONS / SAMPLING_RATE) - CUTOFF) % (int)pow(2, MAX_BTN + 1) == 0) MAX_BTN++;
    char num[(int)((ceil(log10(MAX_BTN)) + 1) * sizeof(char))];
    char errmsg[40] = "BTN must be less than or equal to ";
    sprintf(num, "%d", MAX_BTN);
    strcat(errmsg, num);
    if (BTN > MAX_BTN) {
        error("BTN chosen in config.toml is too large", errmsg);
    }

    // Run simulation using values read from config file then calculate total energy of the last configuration
    double **positions = position_random_generation(N, T, M, FREQUENCY_Z, FREQUENCY_TRANSVERSE);
    double *energies_saved =
        metropolis_hastings(positions, ITERATIONS, N, M, T, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, FREQUENCY_Z,
                            FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT, SAMPLING_RATE);

    printf("Energy of last sampled configuration: %e\n", energies_saved[ITERATIONS / SAMPLING_RATE - 1]);

    double *sliced_energies_saved = slice(energies_saved, CUTOFF, ITERATIONS / SAMPLING_RATE);
    free(energies_saved);
    double *reblocked_energies = reblocking(sliced_energies_saved, (ITERATIONS / SAMPLING_RATE) - CUTOFF, BTN);

    // Calculates the standard errors post equilibration (decided by value of CUTOFF) for different BTN
    double *errors = malloc(MAX_BTN * sizeof(*errors));
    double *mean_energies = reblocking(sliced_energies_saved, (ITERATIONS / SAMPLING_RATE) - CUTOFF, 1);
    free(sliced_energies_saved);
    int size_mean_energies = 0.5 * ((ITERATIONS / SAMPLING_RATE) - CUTOFF);
    errors[0] = calculate_error(mean_energies, size_mean_energies, N);
    for (int i = 1; i < MAX_BTN; i++) {
        double *temp_array = reblocking(mean_energies, size_mean_energies, 1);
        size_mean_energies /= 2;
        mean_energies = realloc(mean_energies, size_mean_energies * sizeof(*mean_energies));
        memcpy(mean_energies, temp_array, size_mean_energies * sizeof(*temp_array));
        free(temp_array);
        errors[i] = calculate_error(mean_energies, size_mean_energies, N);
    }
    free(mean_energies);
    export_1D_array("simulation_energy_data.txt", reblocked_energies,
                    ((ITERATIONS / SAMPLING_RATE) - CUTOFF) / pow(2, BTN));
    export_1D_array("simulation_error_data.txt", errors, MAX_BTN);
    free_2D_array(positions, N);
}

void error(const char *msg, const char *errmsg) {
    fprintf(stderr, "ERROR: %s - %s\n", msg, errmsg);
    exit(-1);
}

void read_config(int *N, int *ITERATIONS, double *T, double *M, double *SIGMA, double *DIPOLE_MOMENT,
                 double *DIPOLE_UNIT_VECTOR, double *FREQUENCY_Z, double *FREQUENCY_TRANSVERSE,
                 double *WALL_REPULSION_COEFFICIENT, int *SAMPLING_RATE, int *BTN, int *CUTOFF) {
    FILE *fp = fopen("config.toml", "r");  // Read and parse toml file
    if (!fp) {
        error("Cannot open config.toml", "");
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
    toml_datum_t data_sampling_rate = toml_int_in(properties, "sampling_rate");
    *SAMPLING_RATE = data_sampling_rate.u.i;
    toml_datum_t blocking_transformation_number = toml_int_in(properties, "blocking_transformation_number");
    *BTN = blocking_transformation_number.u.i;
    toml_datum_t cutoff = toml_int_in(properties, "cutoff");
    *CUTOFF = cutoff.u.i;
    char num[(int)((ceil(log10(*ITERATIONS / *SAMPLING_RATE)) + 1) * sizeof(char))];
    char errmsg[40] = "Cutoff must be less than ";
    sprintf(num, "%d", *ITERATIONS / *SAMPLING_RATE);
    strcat(errmsg, num);
    if (*CUTOFF >= *ITERATIONS / *SAMPLING_RATE) {
        error("Cutoff value chosen in config.toml is too large", errmsg);
    }
    toml_free(conf);  // Free memory
}

double **position_random_generation(int N, double T, double M, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE) {
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);  // Generator type

    // Randomly generate N x 3 array from uniform distribution
    double **array = malloc(N * sizeof(**array));
    const double r_xy = sqrt(6 * BOLTZMANN * T / (M * FREQUENCY_TRANSVERSE * FREQUENCY_TRANSVERSE));
    const double r_z = sqrt(6 * BOLTZMANN * T / (M * FREQUENCY_Z * FREQUENCY_Z));

    for (int i = 0; i < N; i++) {
        array[i] = malloc(3 * sizeof(*array));
        array[i][0] = gsl_ran_flat(r, -r_xy, r_xy);
        array[i][1] = gsl_ran_flat(r, -r_xy, r_xy);
        array[i][2] = gsl_ran_flat(r, -r_z, r_z);
    }
    gsl_rng_free(r);
    return array;
}

double sum(double *a, int D) {
    double result = 0;

    for (int i = 0; i < D; i++) {
        result += a[i];
    }
    return result;
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

// Creates a subset of an array from one given index to another
double *slice(double *a, int from, int until) {
    int size = until - from;
    double *array = malloc(size * sizeof(*array));

    for (int i = 0; i < size; i++) {
        array[i] = a[from + i];
    }
    return array;
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
                        double WALL_REPULSION_COEFFICIENT) {
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
            wall_repulsion += WALL_REPULSION_COEFFICIENT / pow(distance, 12);
        }
    }
    return trapping_potential + (DIPOLE_MOMENT * DIPOLE_MOMENT * dipole_dipole_interaction) + wall_repulsion;
}

// Function that calculates the total potential energy of a given configuration
double calculate_total_energy(double **positions, int N, double M, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR,
                              double FREQUENCY_Z, double FREQUENCY_TRANSVERSE, double WALL_REPULSION_COEFFICIENT) {
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
            total_wall_repulsion += WALL_REPULSION_COEFFICIENT / pow(distance, 12);
        }
    }
    return total_trapping_potential + (DIPOLE_MOMENT * DIPOLE_MOMENT * total_dipole_dipole_interaction) +
           total_wall_repulsion;
}

double *metropolis_hastings(double **positions, int ITERATIONS, int N, double M, double T, double SIGMA,
                            double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z,
                            double FREQUENCY_TRANSVERSE, double WALL_REPULSION_COEFFICIENT, int SAMPLING_RATE) {
    clock_t begin = clock();
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);  // Generator type

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
                                 FREQUENCY_Z, FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT);
            for (int j = 0; j < 3; j++) {
                trial_positions[j] = positions[index][j] + gsl_ran_gaussian_ziggurat(r, SIGMA);
            }
            energy_difference =
                calculate_energy(positions, trial_positions, N, M, index, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR,
                                 FREQUENCY_Z, FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT) -
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
                                           FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT);
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

/* Function that calculates the mean of adjacent energies,
then the mean of those mean energies and so on for a specified number of times */
double *reblocking(double *energies_saved, int size, int BTN) {
    double *array = malloc(size * sizeof(*array));

    memcpy(array, energies_saved, size * sizeof(*energies_saved));
    for (int i = 1; i <= BTN; i++) {
        for (int j = 0; j < size / pow(2, i); j++) {
            array[j] = 0.5 * (array[2 * j] + array[2 * j + 1]);
        }
        array = realloc(array, size / pow(2, i) * sizeof(*array));
    }
    return array;
}

// Function that calculates error in mean energies for each subset of mean energies
double calculate_error(double *mean_energies, int size, int N) {
    double error = 0;
    double mean = sum(mean_energies, size) / size;

    for (int i = 0; i < size; i++) {
        error += (mean_energies[i] - mean) * (mean_energies[i] - mean);
    }
    return sqrt(error / (size - 1)) / N;
}

void export_positions(Double3D *positions_saved) {
    FILE *fp = fopen("simulation_position_data.txt", "w");
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

void export_1D_array(char *file_name, double *array, int size) {
    FILE *fp = fopen(file_name, "w");
    for (int i = 0; i < size; i++) {
        fprintf(fp, "%f\n", array[i]);
    }
    fclose(fp);
    free(array);
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
