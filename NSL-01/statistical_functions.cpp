
/*###################################
File containing the implementation of
different statistical function.
LIST of available functions:
- standard deviation: std_dev();
- block method for uncertainty: prog_average_std_dev_block_method();
- chi squared: chi_sqrd();
###################################*/

#include "statistical_functions.h"

double std_dev(double average, double sqrd_average, int sample_number){
	/*
	Restituisce la deviazione standard normalizzata sulla radice della grandezza del campione
	*/
        if (sample_number == 0){
                return 0;
        } else
                return sqrt( ( sqrd_average - average*average )/sample_number );
};

void prog_average_std_dev_block_method(const string& output_file, double* average_array, double* sqrd_average_array, int number_blocks){
	/*
  Stampa un file nella cartella data/
  contenente la media progressiva e la std-dev delle misure
  sui blocchi
	*/

        ofstream out_file;
        out_file.open(output_file);
        double *prog_average = new double[number_blocks]();         // Definisco il vettore medie progressive
        double *prog_average_sqr = new double[number_blocks]();     // Vettore medie quadre
        double *prog_error = new double[number_blocks]();           // Definisco il vettore con gli errori

	for(int i=0; i < number_blocks; i++){

		for(int j=0; j < i+1; j++){
        		prog_average[i] += average_array[j];
        		prog_average_sqr[i] += sqrd_average_array[j];
		}

		prog_average[i] = prog_average[i]/(i+1);
		prog_average_sqr[i] = prog_average_sqr[i]/(i+1);
		prog_error[i] = std_dev(prog_average[i],prog_average_sqr[i],i);

		out_file << prog_average[i] << " " << prog_error[i] << endl;
	}

	out_file.close();
};

double chi_sqrd(double* observation_vec, double* expected_value_vec, double* variance_vec, int observation_number){
	/*
	Restituisce il chi^2 dei dati
	*/

	double chi_sqrd = 0;
	for(int i=0; i< observation_number; i++){
		chi_sqrd += pow((observation_vec[i] - expected_value_vec[i]),2)/variance_vec[i];
	};

	return chi_sqrd;
};
