#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

double std_dev(double average, double sqrd_average, int sample_number);
void prog_average_std_dev_block_method(const string& name_output_file, double* average_array, double* sqrd_average_array, int number_blocks);
double chi_sqrd(double* observations_vec, double* expected_values_vec, double* variances_vec, int observation_number);
