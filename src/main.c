#include <math.h>
#include <toml.h>
#include <gsl/gsl_randist.h>

typedef struct {
    int m, n, l;
    double* data;
} Double3D;

void read_config(int *N, int *ITERATIONS, double *T, double *M, double *LENGTH, double *SIGMA, double *DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double *FREQUENCY_Z, double *FREQUENCY_TRANSVERSE);
double **position_random_generation(int N, double max);
double sum(double *a, int D);
double magnitude(double *a, int D);
double dot_product(double *a, double *b, int D);
double calculate_energy(double **positions, double *position, int N, double M, int index, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE);
double *calculate_total_energy(double **positions, int N, double M, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE);
void metropolis_hastings(double **positions, int ITERATIONS, int N, double M, int T, double SIGMA, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE);
void export_positions(Double3D *positions_saved);
void export_energies(double *energies_saved, int N);

int main(int argc, char **argv) {
    int N, ITERATIONS;
    double T, M, LENGTH, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR[3], FREQUENCY_Z, FREQUENCY_TRANSVERSE;
    read_config(&N, &ITERATIONS, &T, &M, &LENGTH, &SIGMA, &DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, &FREQUENCY_Z, &FREQUENCY_TRANSVERSE);

    printf("Current variables set in config:\nN: %d\niterations: %d\ntemperature: %f\nmass: %f\nlength: %f\nsigma: %e\ndipole moment magnitude: %f\ndipole unit vector: %f %f %f\nfrequency_z %f\nfrequency_transverse %f\n", N, ITERATIONS, T, M, LENGTH, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR[0], DIPOLE_UNIT_VECTOR[1], DIPOLE_UNIT_VECTOR[2], FREQUENCY_Z, FREQUENCY_TRANSVERSE);

    double **positions = position_random_generation(N, LENGTH);
    metropolis_hastings(positions, ITERATIONS, N, M, T, SIGMA, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, FREQUENCY_Z, FREQUENCY_TRANSVERSE);
    double *energies_saved = calculate_total_energy(positions, N, M, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, FREQUENCY_Z, FREQUENCY_TRANSVERSE);

    printf("Total energy: %e\n", sum(energies_saved, N));
    for (int i = 0; i < N; i++) {
        printf("positions: %f %f %f\n", positions[i][0], positions[i][1], positions[i][2]);
    }

    export_energies(energies_saved, N);

    free(energies_saved);
    for (int i = 0; i < N; i++) {
        free(positions[i]);
    }
    free(positions);
    return 0;
}

void read_config(int *N, int *ITERATIONS, double *T, double *M, double *LENGTH, double *SIGMA, double *DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double *FREQUENCY_Z, double *FREQUENCY_TRANSVERSE) {
    FILE *fp = fopen("config.toml", "r"); // 1. Read and parse toml file
    char errbuf[200];
    toml_table_t *conf = toml_parse_file(fp, errbuf, sizeof(errbuf));
    fclose(fp);
    toml_table_t *properties = toml_table_in(conf, "simulation_properties"); // 2. Traverse to a table.
    toml_datum_t particles = toml_int_in(properties, "particles"); // 3. Extract values
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
    toml_free(conf); // 4. Free memory
}

double **position_random_generation(int N, double max) {
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default); // generator type
    
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

double calculate_energy(double **positions, double *position, int N, double M, int index, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE) {
    const double trapping_potential = 0.5 * M * FREQUENCY_Z * FREQUENCY_Z * ((FREQUENCY_TRANSVERSE / FREQUENCY_Z) * (FREQUENCY_TRANSVERSE / FREQUENCY_Z) * (position[0] * position[0] + position[1] * position[1]) + position[2] * position[2]);
    double dipole_dipole_interaction = 0;
    double displacement[3], distance, vector_term;

    for (int i = 0; i < N; i++) {
        if (i != index) {
            for (int j = 0; j < 3; j++) {
                displacement[j] = position[j] - positions[i][j];
            }
            distance = magnitude(displacement, 3);
            vector_term = dot_product(displacement, DIPOLE_UNIT_VECTOR, 3); 
            dipole_dipole_interaction += 1 / pow(distance, 6) + (DIPOLE_MOMENT * DIPOLE_MOMENT) * (distance * distance - 3 * vector_term * vector_term) / pow(distance, 5);
        }
    }
    // printf("values: %e %f %f\n", M, FREQUENCY_Z, FREQUENCY_TRANSVERSE);
    // printf("positions: %f %f %f\n", position[0], position[1], position[2]);
    // printf("trapping: %e dd int: %e\n", trapping_potential, dipole_dipole_interaction);
    return trapping_potential + dipole_dipole_interaction;
}

double *calculate_total_energy(double **positions, int N, double M, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE) {
    double dipole_dipole_interaction = 0;
    double trapping_potential, displacement[3], distance, vector_term;
    double *energies_saved = malloc(N * sizeof(*energies_saved));

    for (int i = 0; i < N; i++) {
        trapping_potential = 0.5 * M * FREQUENCY_Z * FREQUENCY_Z * ((FREQUENCY_TRANSVERSE / FREQUENCY_Z) * (FREQUENCY_TRANSVERSE / FREQUENCY_Z) * (positions[i][0] * positions[i][0] + positions[i][1] * positions[i][1]) + positions[i][2] * positions[i][2]);
        for (int j = 0; j < N; j++) {
            if (j < i) {
                for (int k = 0; k < 3; k++) {
                    displacement[k] = positions[i][k] - positions[j][k];
                }
                distance = magnitude(displacement, 3);
                vector_term = dot_product(displacement, DIPOLE_UNIT_VECTOR, 3); 
                dipole_dipole_interaction += 1 / pow(distance, 6) + (DIPOLE_MOMENT * DIPOLE_MOMENT) * (distance * distance - 3 * vector_term * vector_term) / pow(distance, 5);
            }
        }
        energies_saved[i] = trapping_potential + dipole_dipole_interaction;
    }
    return energies_saved;
}

void metropolis_hastings(double **positions, int ITERATIONS, int N, double M, int T, double SIGMA, double DIPOLE_MOMENT, double *DIPOLE_UNIT_VECTOR, double FREQUENCY_Z, double FREQUENCY_TRANSVERSE) {
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default); // generator type

    int accepted = 0;
    double trial_positions[3], energy_difference, energy_previous;
    const double kB = 3.167e-6; // Boltzmann constant in Hartree atomic units
    Double3D positions_saved = {ITERATIONS, N, 3};
    positions_saved.data = malloc(positions_saved.m * positions_saved.n * positions_saved.l * sizeof(*positions_saved.data));
    for (int i = 0; i < ITERATIONS; i++) {
        for (int index = 0; index < N; index++) {
            energy_previous = calculate_energy(positions, positions[index], N, M, index, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, FREQUENCY_Z, FREQUENCY_TRANSVERSE);
            for (int j = 0; j < 3; j++) {
                trial_positions[j] = positions[index][j] + gsl_ran_gaussian_ziggurat(r, SIGMA); // make more efficient picking trial configuration for all three directions etc
            }
            energy_difference = calculate_energy(positions, trial_positions, N, M, index, DIPOLE_MOMENT, DIPOLE_UNIT_VECTOR, FREQUENCY_Z, FREQUENCY_TRANSVERSE) - energy_previous;
            if (energy_difference <= 0) {
                accepted++;
                for (int k = 0; k < 3; k++) {
                    positions[index][k] = trial_positions[k];        
                }
            }
            else {
                if (gsl_rng_uniform(r) <= exp(-energy_difference / (kB * T))) {
                    accepted++;
                    for (int k = 0; k < 3; k++) {
                        positions[index][k] = trial_positions[k];        
                    }
                }
            }
            for (int j = 0; j < N; j++) { // change code to only save last 80% of iterations
                for (int k = 0; k < 3; k++) {
                    positions_saved.data[i * (positions_saved.n * positions_saved.l) + j * positions_saved.l + k] = positions[j][k];
                }
            }
        }
    }
    gsl_rng_free(r);
    double percent_accepted = 100 * accepted / (ITERATIONS * N);
    printf("\n\npercent accepted: %f%%\nnumber accepted: %d\n\n\n", percent_accepted, accepted);
    export_positions(&positions_saved);
    free(positions_saved.data);
}

void export_positions(Double3D *positions_saved) {
    FILE *fp = fopen("simulation_position_data.txt", "w");
    for (int i = 0; i < positions_saved->m; i++) {
        for (int j = 0; j < positions_saved->n; j++) {
            fprintf(fp, "%f %f %f\n", positions_saved->data[i * (positions_saved->n * positions_saved->l) + j * positions_saved->l],positions_saved->data[i * (positions_saved->n * positions_saved->l) + j * positions_saved->l + 1], positions_saved->data[i * (positions_saved->n * positions_saved->l) + j * positions_saved->l + 2]);
        }
    }
    fclose(fp);
}

void export_energies(double *energies_saved, int N) {
    FILE *fp = fopen("simulation_energy_data.txt", "w");
    for (int i = 0; i < N; i++) {
        fprintf(fp, "%f\n", energies_saved[i]);
    }
    fclose(fp);
}
