#include <stdio.h>
#include <stdlib.h>

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
