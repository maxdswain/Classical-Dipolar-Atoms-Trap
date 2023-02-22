#pragma once
void error(const char *msg, const char *errmsg);
char *give_file_name(char *file_begin, int count);
char *produce_file_name(char *file_begin);
void export_1D_array(char *file_begin, double *array, int size);
void read_2D_array(char *file_name, double **array, int N);
void free_2D_array(double **array, int size);
