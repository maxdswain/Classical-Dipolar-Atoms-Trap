#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <toml.h>

#include "shared.h"

void read_config(int *ITERATIONS, int *SAMPLING_RATE, int *BTN, int *CUTOFF);
double *read_energies(int ITERATIONS, int SAMPLING_RATE);
double sum(double *a, int D);
double *slice(double *a, int from, int until);
double *reblocking(double *energies, int size, int BTN);
double calculate_error(double *mean_energies, int size);

int main(int argc, char **argv) {
    int ITERATIONS, SAMPLING_RATE, BTN, CUTOFF;
    read_config(&ITERATIONS, &SAMPLING_RATE, &BTN, &CUTOFF);
    double *energies = read_energies(ITERATIONS, SAMPLING_RATE);
    double *sliced_energies = slice(energies, CUTOFF, ITERATIONS / SAMPLING_RATE);

    int MAX_BTN = 2;  // ((ITERATIONS / SAMPLING_RATE) - CUTOFF) / 2^MAX_BTN must be an integer
    while (((ITERATIONS / SAMPLING_RATE) - CUTOFF) % (int)pow(2, MAX_BTN + 1) == 0) MAX_BTN++;
    char num[(int)((ceil(log10(MAX_BTN)) + 1) * sizeof(char))];
    char errmsg[40] = "BTN must be less than or equal to ";
    sprintf(num, "%d", MAX_BTN);
    strcat(errmsg, num);
    if (BTN > MAX_BTN) {
        error("BTN chosen in input.toml is too large", errmsg);
    }

    // Calculates the standard errors post equilibration (decided by value of CUTOFF) for different BTN
    double *errors = malloc(MAX_BTN * sizeof(*errors));
    double *mean_energies = reblocking(sliced_energies, (ITERATIONS / SAMPLING_RATE) - CUTOFF, 1);
    int size_mean_energies = 0.5 * ((ITERATIONS / SAMPLING_RATE) - CUTOFF);
    errors[0] = calculate_error(mean_energies, size_mean_energies);
    for (int i = 1; i < MAX_BTN; i++) {
        double *temp_array = reblocking(mean_energies, size_mean_energies, 1);
        size_mean_energies /= 2;
        mean_energies = realloc(mean_energies, size_mean_energies * sizeof(*mean_energies));
        memcpy(mean_energies, temp_array, size_mean_energies * sizeof(*temp_array));
        free(temp_array);
        errors[i] = calculate_error(mean_energies, size_mean_energies);
    }
    free(energies);
    free(sliced_energies);
    export_1D_array("mean_energies", mean_energies, size_mean_energies);
    export_1D_array("error_data", errors, MAX_BTN);
}

void read_config(int *ITERATIONS, int *SAMPLING_RATE, int *BTN, int *CUTOFF) {
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
    toml_datum_t repetitions = toml_int_in(properties, "repetitions");        // Extract values from table
    *ITERATIONS = repetitions.u.i;
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
        error("Cutoff value chosen in input.toml is too large", errmsg);
    }
    toml_free(conf);  // Free memory
}

double *read_energies(int ITERATIONS, int SAMPLING_RATE) {
    FILE *fp = fopen("energy_data.out", "r");
    if (!fp) {
        error("Cannot find energy_data.out", strerror(errno));
    }
    double *energies = malloc(ITERATIONS / SAMPLING_RATE * sizeof(*energies));
    for (int i = 0; i < ITERATIONS / SAMPLING_RATE; i++) {
        if (fscanf(fp, "%lf", &energies[i]) != 1) {
            char num[(int)((ceil(log10(ITERATIONS / SAMPLING_RATE)) + 1) * sizeof(char))];
            char errmsg[100] = "There should be a double per line for ";
            sprintf(num, "%d", ITERATIONS / SAMPLING_RATE);
            strcat(errmsg, num);
            strcat(errmsg, " lines");
            error("energy_data.out is not in the correct format", errmsg);
        }
    }
    fclose(fp);
    return energies;
}

double sum(double *a, int D) {
    double result = 0;

    for (int i = 0; i < D; i++) {
        result += a[i];
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

/* Function that calculates the mean of adjacent energies,
then the mean of those mean energies and so on for a specified number of times */
double *reblocking(double *energies, int size, int BTN) {
    double *array = malloc(size * sizeof(*array));

    memcpy(array, energies, size * sizeof(*energies));
    for (int i = 1; i <= BTN; i++) {
        for (int j = 0; j < size / pow(2, i); j++) {
            array[j] = 0.5 * (array[2 * j] + array[2 * j + 1]);
        }
        array = realloc(array, size / pow(2, i) * sizeof(*array));
    }
    return array;
}

// Function that calculates error in mean energies for each subset of mean energies
double calculate_error(double *mean_energies, int size) {
    double error = 0;
    double mean = sum(mean_energies, size) / size;

    for (int i = 0; i < size; i++) {
        error += (mean_energies[i] - mean) * (mean_energies[i] - mean);
    }
    return sqrt(error / size);
}
