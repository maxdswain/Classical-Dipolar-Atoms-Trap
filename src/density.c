#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <toml.h>

#include "shared.h"

void read_config(int *N, int *BINS_X, int *BINS_Y, int *BINS_Z);
double max_position(double **configuration, int index, int size);
double min_position(double **configuration, int index, int size);
double *calculate_density(double **configuration, int N, int BINS_X, int BINS_Y, int BINS_Z);
double *calculate_pair_density(double **configuration, int N, int BINS_X, int BINS_Y, int BINS_Z);

int main(int argc, char **argv) {
    // Initialise and read needed values from input file.
    int N, BINS_X, BINS_Y, BINS_Z;
    read_config(&N, &BINS_X, &BINS_Y, &BINS_Z);

    // Extract last iteration from positions output file.
    char num[(int)((ceil(log10(N)) + 1) * sizeof(char))];
    char command[50] = "tail -n ";
    sprintf(num, "%d", N);
    strcat(command, num);
    strcat(command, " position_data.out > configuration.out");
    if (system(command)) {
        error("Cannot run tail system command", strerror(errno));
    }
    double **configuration = malloc(N * sizeof(**configuration));
    for (int i = 0; i < N; i++) {
        configuration[i] = malloc(3 * sizeof(*configuration));
    }
    read_2D_array("configuration.out", configuration, N);

    double *density = calculate_density(configuration, N, BINS_X, BINS_Y, BINS_Z);
    export_1D_array("density.out", density, BINS_X * BINS_Y * BINS_Z);

    double *pair_density = calculate_pair_density(configuration, N, BINS_X, BINS_Y, BINS_Z);
    export_1D_array("pair_density.out", pair_density, BINS_X * BINS_Y * BINS_Z * BINS_X * BINS_Y * BINS_Z);

    free_2D_array(configuration, N);
}

void read_config(int *N, int *BINS_X, int *BINS_Y, int *BINS_Z) {
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
    toml_datum_t x_bins = toml_int_in(properties, "bins_x");
    *BINS_X = x_bins.u.i;
    toml_datum_t y_bins = toml_int_in(properties, "bins_y");
    *BINS_Y = y_bins.u.i;
    toml_datum_t z_bins = toml_int_in(properties, "bins_z");
    *BINS_Z = z_bins.u.i;
    toml_free(conf);  // Free memory
}

double max_position(double **configuration, int index, int size) {
    double max = configuration[0][index];
    for (int i = 0; i < size; i++) {
        if (max < configuration[i][index]) {
            max = configuration[i][index];
        }
    }
    return max;
}

double min_position(double **configuration, int index, int size) {
    double min = configuration[0][index];
    for (int i = 0; i < size; i++) {
        if (min > configuration[i][index]) {
            min = configuration[i][index];
        }
    }
    return min;
}

double *calculate_density(double **configuration, int N, int BINS_X, int BINS_Y, int BINS_Z) {
    // Calculate maximum-minimum x y z positions - use an area a little bit larger than min to max
    double min_x = min_position(configuration, 0, N);
    double bin_length_x = 1.1 * (max_position(configuration, 0, N) - min_x) / (double)BINS_X;
    double min_y = min_position(configuration, 1, N);
    double bin_length_y = 1.1 * (max_position(configuration, 1, N) - min_y) / (double)BINS_Y;
    double min_z = min_position(configuration, 2, N);
    double bin_length_z = 1.1 * (max_position(configuration, 2, N) - min_z) / (double)BINS_Z;
    // Loop over all bin widths and calculate density
    double *density = malloc(BINS_X * BINS_Y * BINS_Z * sizeof(*density));
    for (int i = 0; i < BINS_X; i++) {
        for (int j = 0; j < BINS_Y; j++) {
            for (int k = 0; k < BINS_Z; k++) {
                int counter = 0;
                for (int n = 0; n < N; n++) {
                    // If atom in bin add one to counter
                    double x = configuration[n][0];
                    double y = configuration[n][1];
                    double z = configuration[n][2];
                    if ((x >= 1.05 * min_x + i * bin_length_x && x < 1.05 * min_x + (i + 1) * bin_length_x) &&
                        (y >= 1.05 * min_y + j * bin_length_y && y < 1.05 * min_y + (j + 1) * bin_length_y) &&
                        (z >= 1.05 * min_z + k * bin_length_z && z < 1.05 * min_z + (k + 1) * bin_length_z)) {
                        counter++;
                    }
                }
                density[i * (BINS_Y * BINS_Z) + j * BINS_Z + k] = counter;
            }
        }
    }
    return density;
}

double *calculate_pair_density(double **configuration, int N, int BINS_X, int BINS_Y, int BINS_Z) {
    // Calculate maximum-minimum x y z positions - use an area a little bit larger than min to max
    double min_x = min_position(configuration, 0, N);
    double bin_length_x = 1.1 * (max_position(configuration, 0, N) - min_x) / (double)BINS_X;
    double min_y = min_position(configuration, 1, N);
    double bin_length_y = 1.1 * (max_position(configuration, 1, N) - min_y) / (double)BINS_Y;
    double min_z = min_position(configuration, 2, N);
    double bin_length_z = 1.1 * (max_position(configuration, 2, N) - min_z) / (double)BINS_Z;
    // Loop over all bin widths and calculate pair density
    double *pair_density = malloc(BINS_X * BINS_Y * BINS_Z * BINS_X * BINS_Y * BINS_Z * sizeof(*pair_density));
    for (int i = 0; i < BINS_X; i++) {
        for (int j = 0; j < BINS_Y; j++) {
            for (int k = 0; k < BINS_Z; k++) {
                for (int l = 0; l < BINS_X; l++) {
                    for (int m = 0; m < BINS_Y; m++) {
                        for (int n = 0; n < BINS_Z; n++) {
                            int counter = 0;
                            // Pairwise sum over both bins, then counting pairs
                            for (int a = 0; a < N; a++) {
                                double x = configuration[a][0];
                                double y = configuration[a][1];
                                double z = configuration[a][2];
                                if ((x >= 1.05 * min_x + i * bin_length_x &&
                                     x < 1.05 * min_x + (i + 1) * bin_length_x) &&
                                    (y >= 1.05 * min_y + j * bin_length_y &&
                                     y < 1.05 * min_y + (j + 1) * bin_length_y) &&
                                    (z >= 1.05 * min_z + k * bin_length_z &&
                                     z < 1.05 * min_z + (k + 1) * bin_length_z)) {
                                    for (int b = 0; b < a; b++) {
                                        x = configuration[b][0];
                                        y = configuration[b][1];
                                        z = configuration[b][2];
                                        if ((x >= 1.05 * min_x + l * bin_length_x &&
                                             x < 1.05 * min_x + (l + 1) * bin_length_x) &&
                                            (y >= 1.05 * min_y + m * bin_length_y &&
                                             y < 1.05 * min_y + (m + 1) * bin_length_y) &&
                                            (z >= 1.05 * min_z + n * bin_length_z &&
                                             z < 1.05 * min_z + (n + 1) * bin_length_z)) {
                                            counter++;
                                        }
                                    }
                                }
                            }
                            // Store pair density by storing counts for each i, j, k, l, m, n in 6D array in 1D
                            pair_density[i * (BINS_Y * BINS_Z * BINS_X * BINS_Y * BINS_Z) +
                                         j * (BINS_Z * BINS_X * BINS_Y * BINS_Z) + k * (BINS_X * BINS_Y * BINS_Z) +
                                         l * (BINS_Y * BINS_Z) + m * BINS_Z + n] = counter;
                        }
                    }
                }
            }
        }
    }
    return pair_density;
}
