#include <math.h>
#include <toml.h>
#include <gsl/gsl_randist.h>

void read_config(int* N, int* iterations, int* T, double* m, double* length, double* sigma, double* dipole_moment, double *frequency_z, double *frequency_transverse);
double **position_random_generation(int N, double max);
double magnitude(const double *a, int D);
double dot_product(const double *a, const double *b, int D);
double calculate_energy(double **positions, double *position, int N, double m, int index, double dipole_moment, double frequency_z, double frequency_transverse);
double calculate_total_energy(double **positions, int N, double m, double dipole_moment, double frequency_z, double frequency_transverse);
void particle_metropolis_hastings(double **positions, int iterations, int N, double m, int T, double sigma, double dipole_moment, double frequency_z, double frequency_transverse);
void export_positions(double ***positions_saved, int iterations, int N);

int main(int argc, char **argv) {
    int N, iterations, T;
    double m, length, sigma, dipole_moment, frequency_z, frequency_transverse;
    read_config(&N, &iterations, &T, &m, &length, &sigma, &dipole_moment, &frequency_z, &frequency_transverse);

    printf("Current variables set in config:\nN: %d\niterations: %d\ntemperature: %d\nmass: %e\nlength: %f\nsigma: %e\nmagnetic dipole moment: %e\nfrequency_z %f\nfrequency_transverse %f\n", N, iterations, T, m, length, sigma, dipole_moment, frequency_z, frequency_transverse);

    double **positions = position_random_generation(N, length);
    particle_metropolis_hastings(positions, iterations, N, m, T, sigma, dipole_moment, frequency_z, frequency_transverse);

    printf("Total energy: %e\n", calculate_total_energy(positions, N, m, dipole_moment, frequency_z, frequency_transverse));

    for (int i = 0; i < N; i++) { // Tests
        printf("positions: %f %f %f\n", positions[i][0], positions[i][1], positions[i][2]);
    }
    for (int i = 0; i < N; i++) {
        free(positions[i]);
    }
    free(positions);
    return 0;
}

void read_config(int* N, int* iterations, int* T, double* m, double* length, double* sigma, double* dipole_moment, double *frequency_z, double *frequency_transverse) {
    FILE *fp = fopen("config.toml", "r"); // 1. Read and parse toml file
    char errbuf[200];
    toml_table_t *conf = toml_parse_file(fp, errbuf, sizeof(errbuf));
    fclose(fp);
    toml_table_t *properties = toml_table_in(conf, "simulation_properties"); // 2. Traverse to a table.
    toml_datum_t particles = toml_int_in(properties, "particles"); // 3. Extract values
    *N = particles.u.i;
    toml_datum_t repetitions = toml_int_in(properties, "repetitions");
    *iterations = repetitions.u.i;
    toml_datum_t temperature = toml_int_in(properties, "temperature");
    *T = temperature.u.i;
    toml_datum_t mass = toml_double_in(properties, "mass");
    *m = mass.u.d;
    toml_datum_t box_length = toml_double_in(properties, "box_length");
    *length = box_length.u.d;
    toml_datum_t trial_sigma = toml_double_in(properties, "sigma");
    *sigma = trial_sigma.u.d;
    toml_datum_t magnetic_dipole_moment = toml_double_in(properties, "magnetic_dipole_moment");
    *dipole_moment = magnetic_dipole_moment.u.d;
    toml_datum_t trapping_frequency_z = toml_double_in(properties, "trapping_frequency_z");
    *frequency_z = trapping_frequency_z.u.d;
    toml_datum_t trapping_frequency_transverse = toml_double_in(properties, "trapping_frequency_transverse");
    *frequency_transverse = trapping_frequency_transverse.u.d;
    toml_free(conf); // 4. Free memory
}

double **position_random_generation(int N, double max) {
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default); // generator type
    
    double **array = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        array[i] = (double *)malloc(3 * sizeof(double));
        for (int j = 0; j < 3; j++) {
            array[i][j] = gsl_rng_uniform_int(r, max);
        }
    }
    gsl_rng_free(r);
    return array;
}

double magnitude(const double *a, int D) {
    double result = 0.0;
    for (int i = 0; i < D; i++) {
        result += a[i] * a[i];
    }
    return sqrt(result);
}

double dot_product(const double *a, const double *b, int D) {
    double result = 0.0;
    for (int i = 0; i < D; i++) {
        result += a[i] * b[i];
    }
    return result;
}

double calculate_energy(double **positions, double *position, int N, double m, int index, double dipole_moment, double frequency_z, double frequency_transverse) {
    const double mu_zero = 1.25663706212e-6; // vacuum permeability in SI units
    const double trapping_potential = 0.5 * m * frequency_z * frequency_z * ((frequency_transverse / frequency_z) * (frequency_transverse / frequency_z) * (position[0] * position[0] + position[1] * position[1]) + position[2] * position[2]);
    const double dipole_unit_vector[3] = {8 / sqrt(74), 3 / sqrt(74), 1 / sqrt(74)};
    double dipole_dipole_interaction = 0;
    double displacement[3], distance, vector_term;

    for (int i = 0; i < N; i++) {
        if (i != index) {
            for (int j = 0; j < 3; j++) {
                displacement[j] = position[j] - positions[i][j];
            }
            distance = magnitude(displacement, 3);
            vector_term = dot_product(displacement, dipole_unit_vector, 3); 
            dipole_dipole_interaction += (mu_zero * dipole_moment * dipole_moment) * (distance * distance - 3 * vector_term * vector_term) / (4 * M_PI * pow(distance, 5));
        }
    }
    // printf("values: %e %f %f\n", m, frequency_z, frequency_transverse);
    // printf("positions: %f %f %f\n", position[0], position[1], position[2]);
    // printf("trapping: %e dd int: %e\n", trapping_potential, dipole_dipole_interaction);
    return trapping_potential + dipole_dipole_interaction;
}

double calculate_total_energy(double **positions, int N, double m, double dipole_moment, double frequency_z, double frequency_transverse) {
    const double mu_zero = 1.25663706212e-6; // vacuum permeability in SI units
    const double dipole_unit_vector[3] = {8 / sqrt(74), 3 / sqrt(74), 1 / sqrt(74)};
    double total_dipole_dipole_interaction = 0;
    double displacement[3], distance, vector_term, total_trapping_potential;

    for (int i = 0; i < N; i++) {
        total_trapping_potential = 0.5 * m * frequency_z * frequency_z * ((frequency_transverse / frequency_z) * (frequency_transverse / frequency_z) * (positions[i][0] * positions[i][0] + positions[i][1] * positions[i][1]) + positions[i][2] * positions[i][2]);
        for (int j = 0; j < N; j++) {
            if (j < i) {
                for (int k = 0; k < 3; k++) {
                    displacement[k] = positions[i][k] - positions[j][k];
                }
                distance = magnitude(displacement, 3);
                vector_term = dot_product(displacement, dipole_unit_vector, 3); 
                total_dipole_dipole_interaction += (mu_zero * dipole_moment * dipole_moment) * (distance * distance - 3 * vector_term * vector_term) / (4 * M_PI * pow(distance, 5));
            }
        }
    }
    return 0.5 * (total_trapping_potential + total_dipole_dipole_interaction);
}

void particle_metropolis_hastings(double **positions, int iterations, int N, double m, int T, double sigma, double dipole_moment, double frequency_z, double frequency_transverse) {
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default); // generator type

    int accepted = 0;
    double trial_positions[3], energy_difference, energy_previous;
    const double kB = 1.380649e-23; // Boltzmann constant in SI units
    double ***positions_saved = (double ***)malloc(iterations * sizeof(double **));
    for (int i = 0; i < iterations; i++) {
        for (int index = 0; index < N; index++) {
            energy_previous = calculate_energy(positions, positions[index], N, m, index, dipole_moment, frequency_z, frequency_transverse);
            for (int j = 0; j < 3; j++) {
                trial_positions[j] = positions[index][j] + gsl_ran_gaussian_ziggurat(r, sigma); // make more efficient picking trial configuration for all three directions etc
            }
            energy_difference = calculate_energy(positions, trial_positions, N, m, index, dipole_moment, frequency_z, frequency_transverse) - energy_previous;
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
            positions_saved[i] = (double **)malloc(N * sizeof(double *)); // change code to only save last 80% of iterations
            for (int j = 0; j < N; j++) {
                positions_saved[i][j] = (double *)malloc(3 * sizeof(double));
                for (int k = 0; k < 3; k++) {
                    positions_saved[i][j][k] = positions[j][k];
                }
            }
        }
    }
    double percent_accepted = 100 * accepted / (iterations * N);
    printf("\n\npercent accepted: %f\nnumber accepted: %d\n\n\n", percent_accepted, accepted);
    export_positions(positions_saved, iterations, N);
    for (int i = 0; i < iterations; i++) {
        for (int j = 0; j < N; j++) {
            free(positions_saved[i][j]);
        }
        free(positions_saved[i]);
    }
    free(positions_saved);
    gsl_rng_free(r);
}

void export_positions(double ***positions_saved, int iterations, int N) {
    FILE *fp = fopen("simulation_data.txt", "w");
    for (int i = 0; i < iterations; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(fp, "%f %f %f\n", positions_saved[i][j][0], positions_saved[i][j][1], positions_saved[i][j][2]);
        }
    }
    fclose(fp);
}
