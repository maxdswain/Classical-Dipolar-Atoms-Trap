#include <math.h>
#include <time.h>

#include <gsl/gsl_randist.h>
#include <toml.h>

typedef struct {
    int m, n, l;
    double *data;
} Double3D;

void read_config(int *N, int *ITERATIONS, double *T, double *M, double *LENGTH, double *SIGMA, double *DIPOLE_MOMENT,
                 double *DIPOLE_UNIT_VECTOR, double *FREQUENCY_Z, double *FREQUENCY_TRANSVERSE,
                 double *WALL_REPULSION_COEFFICIENT, int *SAMPLING_RATE, int *BTN);
double **position_random_generation(int N, double max);
double sum(double *a, int D);
double magnitude(double *a, int D);
double dot_product(double *a, double *b, int D);
double calculate_energy(double **positions, double *position, int N, double M, int index, double DIPOLE_MOMENT,
                        double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE,
                        double WALL_REPULSION_COEFFICIENT);
double calculate_total_energy(double **positions, int N, double M, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR,
                              double FREQUENCY_Z, double FREQUENCY_TRANSVERSE, double WALL_REPULSION_COEFFICIENT);
void metropolis_hastings(double **positions, int ITERATIONS, int N, double M, int T, double SIGMA, double DIPOLE_MOMENT,
                         double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE,
                         double WALL_REPULSION_COEFFICIENT, int SAMPLING_RATE, int BTN);
double *reblocking(double *energies_saved, int size, int BTN);
void export_positions(Double3D *positions_saved);
void export_energies(double *mean_energies, int size);
void progress_bar(double progress, double time);

int main(int argc, char **argv) {
    int N, ITERATIONS, SAMPLING_RATE, BTN;
    double T, M, LENGTH, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR[3], FREQUENCY_Z, FREQUENCY_TRANSVERSE,
        WALL_REPULSION_COEFFICIENT;
    read_config(&N, &ITERATIONS, &T, &M, &LENGTH, &SIGMA, &DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, &FREQUENCY_Z,
                &FREQUENCY_TRANSVERSE, &WALL_REPULSION_COEFFICIENT, &SAMPLING_RATE, &BTN);

    printf(
        "Current variables set in config:\nN: %d\niterations: %d\ntemperature: %f\nmass: %f\nlength: %f\nsigma: "
        "%e\ndipole moment magnitude: %f\ndipole unit vector: %f %f %f\nfrequency_z %f\nfrequency_transverse %f\nhard "
        "wall repulsion coefficient %f\n\n",
        N, ITERATIONS, T, M, LENGTH, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR[0], DIPOLE_UNIT_VECTOR[1],
        DIPOLE_UNIT_VECTOR[2], FREQUENCY_Z, FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT);

    double **positions = position_random_generation(N, LENGTH);
    metropolis_hastings(positions, ITERATIONS, N, M, T, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, FREQUENCY_Z,
                        FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT, SAMPLING_RATE, BTN);
    double total_energy = calculate_total_energy(positions, N, M, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, FREQUENCY_Z,
                                                 FREQUENCY_TRANSVERSE, WALL_REPULSION_COEFFICIENT);

    printf("Total energy: %e\n", total_energy);

    for (int i = 0; i < N; i++) {
        free(positions[i]);
    }
    free(positions);
    return 0;
}

void read_config(int *N, int *ITERATIONS, double *T, double *M, double *LENGTH, double *SIGMA, double *DIPOLE_MOMENT,
                 double *DIPOLE_UNIT_VECTOR, double *FREQUENCY_Z, double *FREQUENCY_TRANSVERSE,
                 double *WALL_REPULSION_COEFFICIENT, int *SAMPLING_RATE, int *BTN) {
    FILE *fp = fopen("config.toml", "r");  // 1. Read and parse toml file
    char errbuf[200];
    toml_table_t *conf = toml_parse_file(fp, errbuf, sizeof(errbuf));
    fclose(fp);
    toml_table_t *properties = toml_table_in(conf, "simulation_properties");  // 2. Traverse to a table.
    toml_datum_t particles = toml_int_in(properties, "particles");            // 3. Extract values
    *N = particles.u.i;
    toml_datum_t repetitions = toml_int_in(properties, "repetitions");
    *ITERATIONS = repetitions.u.i;
    toml_datum_t temperature = toml_double_in(properties, "temperature");
    *T = temperature.u.d;
    toml_datum_t mass = toml_double_in(properties, "mass");
    *M = mass.u.d;
    toml_datum_t box_length = toml_double_in(properties, "box_length");
    *LENGTH = box_length.u.d;
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
    toml_datum_t hard_wall_repulsion_coefficient = toml_double_in(properties, "hard_wall_repulsion_coefficient");
    *WALL_REPULSION_COEFFICIENT = hard_wall_repulsion_coefficient.u.d;
    toml_datum_t data_sampling_rate = toml_int_in(properties, "sampling_rate");
    *SAMPLING_RATE = data_sampling_rate.u.i;
    toml_datum_t blocking_transformation_number = toml_int_in(properties, "blocking_transformation_number");
    *BTN = blocking_transformation_number.u.i;
    toml_free(conf);  // 4. Free memory
}

double **position_random_generation(int N, double max) {
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);  // generator type

    double **array = malloc(N * sizeof(*array));
    for (int i = 0; i < N; i++) {
        array[i] = malloc(3 * sizeof(array[0]));
        for (int j = 0; j < 3; j++) {
            array[i][j] = gsl_rng_uniform_int(r, max);
        }
    }
    gsl_rng_free(r);
    return array;
}

double sum(double *a, int D) {
    double result = 0.0;
    for (int i = 0; i < D; i++) {
        result += a[i];
    }
    return result;
}

double magnitude(double *a, int D) {
    double result = 0.0;
    for (int i = 0; i < D; i++) {
        result += a[i] * a[i];
    }
    return sqrt(result);
}

double dot_product(double *a, double *b, int D) {
    double result = 0.0;
    for (int i = 0; i < D; i++) {
        result += a[i] * b[i];
    }
    return result;
}

double calculate_energy(double **positions, double *position, int N, double M, int index, double DIPOLE_MOMENT,
                        double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE,
                        double WALL_REPULSION_COEFFICIENT) {
    const double trapping_potential =
        0.5 * M *
        (FREQUENCY_TRANSVERSE * FREQUENCY_TRANSVERSE * (position[0] * position[0] + position[1] * position[1]) +
         FREQUENCY_Z * FREQUENCY_Z * position[2] * position[2]);
    double dipole_dipole_interaction = 0;
    double hard_wall_repulsion = 0;
    double displacement[3], distance, vector_term;

    for (int i = 0; i < N; i++) {
        if (i != index) {
            for (int j = 0; j < 3; j++) {
                displacement[j] = position[j] - positions[i][j];
            }
            distance = magnitude(displacement, 3);
            vector_term = dot_product(displacement, DIPOLE_UNIT_VECTOR, 3);
            dipole_dipole_interaction += (distance * distance - 3 * vector_term * vector_term) / pow(distance, 5);
            hard_wall_repulsion += WALL_REPULSION_COEFFICIENT / pow(distance, 6);
        }
    }
    return trapping_potential + (DIPOLE_MOMENT * DIPOLE_MOMENT * dipole_dipole_interaction) + hard_wall_repulsion;
}

double calculate_total_energy(double **positions, int N, double M, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR,
                              double FREQUENCY_Z, double FREQUENCY_TRANSVERSE, double WALL_REPULSION_COEFFICIENT) {
    double total_dipole_dipole_interaction = 0;
    double total_trapping_potential = 0;
    double total_hard_wall_repulsion = 0;
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
            total_hard_wall_repulsion += WALL_REPULSION_COEFFICIENT / pow(distance, 6);
        }
    }
    return total_trapping_potential + (DIPOLE_MOMENT * DIPOLE_MOMENT * total_dipole_dipole_interaction) +
           total_hard_wall_repulsion;
}

void metropolis_hastings(double **positions, int ITERATIONS, int N, double M, int T, double SIGMA, double DIPOLE_MOMENT,
                         double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE,
                         double WALL_REPULSION_COEFFICIENT, int SAMPLING_RATE, int BTN) {
    clock_t begin = clock();

    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);  // generator type

    int accepted = 0;
    double trial_positions[3], energy_difference, energy_previous;
    const double kB = 3.167e-6;  // Boltzmann constant in Hartree atomic units
    double *energies_saved = malloc(ITERATIONS / SAMPLING_RATE * sizeof(*energies_saved));
    Double3D positions_saved = {ITERATIONS / SAMPLING_RATE, N, 3};
    positions_saved.data =
        malloc(positions_saved.m * positions_saved.n * positions_saved.l * sizeof(*positions_saved.data));
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
                for (int k = 0; k < 3; k++) {
                    positions[index][k] = trial_positions[k];
                }
            } else if (gsl_rng_uniform(r) <= exp(-energy_difference / (kB * T))) {
                accepted++;
                for (int k = 0; k < 3; k++) {
                    positions[index][k] = trial_positions[k];
                }
            }
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
        clock_t end = clock();
        double time = (double)(end - begin) / CLOCKS_PER_SEC;
        progress_bar(((double)i + 1) / (double)ITERATIONS, time);
    }
    gsl_rng_free(r);
    double percent_accepted = 100 * accepted / (ITERATIONS * N);
    printf("\n\npercent accepted: %.2f%%\nnumber accepted: %d\n\n\n", percent_accepted, accepted);
    double *mean_energies = reblocking(energies_saved, ITERATIONS / SAMPLING_RATE, BTN);
    free(energies_saved);
    export_energies(mean_energies, (ITERATIONS / SAMPLING_RATE) / pow(2, BTN));
    free(mean_energies);
    export_positions(&positions_saved);
    free(positions_saved.data);
}

double *reblocking(double *energies_saved, int size, int BTN) {
    double *array = malloc(size * sizeof(*array));
    double *temp_array = malloc(0.5 * size * sizeof(*temp_array));
    for (int i = 0; i < size; i++) {
        array[i] = energies_saved[i];
    }
    for (int i = 1; i <= BTN; i++) {
        for (int j = 0; j < size / pow(2, i); j++) {
            temp_array[j] = 0.5 * (array[2 * j] + array[2 * j + 1]);
        }
        for (int j = 0; j < size / pow(2, i); j++) {
            array[j] = temp_array[j];
        }
    }
    free(temp_array);
    return array;
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
}

void export_energies(double *mean_energies, int size) {
    FILE *fp = fopen("simulation_energy_data.txt", "w");
    for (int i = 0; i < size; i++) {
        fprintf(fp, "%f\n", mean_energies[i]);
    }
    fclose(fp);
}

void progress_bar(double progress, double time) {
    int barWidth = 70;

    printf("\r|");
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; i++) {
        if (i < pos) {
            printf("█");
        } else if (i == pos) {
            printf(">");
        } else {
            printf("-");
        }
    }
    if (progress == 1.0) {
        printf("| %.2f %% Total time taken: %.2fs", progress * 100.0, time);
        fflush(stdout);
    } else {
        printf("| %.2f %%", progress * 100.0);
        fflush(stdout);
    }
}
