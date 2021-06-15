#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"

using namespace std;

int main(int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;

	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;


	int M = 1E8;			// Total number of throws
	int N = 1E2;                 	// Number of blocks
	int L = static_cast<int>(M/N);  // Number of throws in each block, please use for M a multiple of N
	double *random_vec = new double[M]();		// Define random vector

	for(int i=0; i < M; i++){	// Load the vector with random number distributed uniformly
		random_vec[i] = rnd.Rannyu();
	}

	// Estimate of <r>
	double *average1 = new double[N]();		// Define average vector
	double *average_sqr1 = new double[N]();		// Define average squared vector
	for(int i=0; i < N; i++){       // Compute the average of my observable and the aveˆ2 to calculate the variance
		double sum = 0;
		for(int j=0; j < L; j++){
			int k = j + i*L;
			sum = sum + random_vec[k];
		}
		average1[i] = sum/L;
		average_sqr1[i] = pow(average1[i],2);
	}
	prog_average_std_dev_block_method("results/EX011(1).dat", average1, average_sqr1, N);
	// Estimate of <(sigma_r)^2>
        double *average2 = new double[N]();             // Define average vector
        double *average_sqr2 = new double[N]();         // Define average squared vector
	for(int i=0; i < N; i++){       // Compute the average of my observable and the aveˆ2 to calculate the variance
                 double sum = 0;
                 for(int j=0; j < L; j++){
                         int k = j + i*L;
                         sum = sum + pow((random_vec[k]-0.5),2);
                 }
                 average2[i] = sum/L;
                 average_sqr2[i] = pow(average2[i],2);
         }
         prog_average_std_dev_block_method("results/EX011(2).dat", average2, average_sqr2, N);

	// Estimate chi^2 dividing [0,1] in sub-intervals intervals=10^2 part
	int intervals = 100;
	double interval_prob = 1./intervals;		// Probability to sit in a specific sub-interval due to uniform distr
	int throws = static_cast<int>(L/intervals); 	// Throws in each sub-interval

	double expected_value = throws/intervals;	// Expected value due to uniform distribution
	double *expected_value_vec = new double[intervals]();
	for(int i=0; i< intervals; i++) expected_value_vec[i] = expected_value;

	double *observations = new double[intervals]();

        ofstream out_file;
        out_file.open("results/chi_sqrd.dat");
	for(int t=0; t < N; t++){

		for(int i=0; i < intervals; i++){
			int n_hits = 0;
			for(int j=0; j < throws ;j++){
				int k = j + i*throws + t*throws*N;
				if( random_vec[k] >= (0 + i*interval_prob) && random_vec[k] <= (i*interval_prob + interval_prob)){
					n_hits += 1;
				}
			observations[i] = n_hits;
			}
		};
    double chi;
		chi = chi_sqrd(observations, expected_value_vec, expected_value_vec, intervals);
                out_file << chi << endl;
	};
        out_file.close();

	return 0;
}
