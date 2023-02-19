#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void error(const char *msg, const char *errmsg) {
    fprintf(stderr, "ERROR: %s - %s\n", msg, errmsg);
    exit(-1);
}

void export_1D_array(char *file_name, double *array, int size) {
    FILE *fp = fopen(file_name, "w");
    for (int i = 0; i < size; i++) {
        fprintf(fp, "%f\n", array[i]);
    }
    fclose(fp);
    free(array);
}

void read_2D_array(char *file_name, double **array, int N) {
    FILE *fp = fopen(file_name, "r");
    if (!fp) {
        char errmsg[40] = "Cannot find ";
        strcat(errmsg, file_name);
        error(errmsg, strerror(errno));
    }
    for (int i = 0; i < N; i++) {
        if (fscanf(fp, "%lf %lf %lf", &array[i][0], &array[i][1], &array[i][2]) != 3) {
            char errmsg[100];
            strcpy(errmsg, file_name);
            strcat(errmsg, " is not in the correct format");
            char num[(int)((ceil(log10(N)) + 1) * sizeof(char))];
            char errmsg2[100] = "There should be 3 doubles per line for ";
            sprintf(num, "%d", N);
            strcat(errmsg2, num);
            strcat(errmsg2, " lines");
            error(errmsg, errmsg2);
        }
    }
    fclose(fp);
}

void free_2D_array(double **array, int size) {
    for (int i = 0; i < size; i++) {
        free(array[i]);
    }
    free(array);
}
