#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"

using namespace std;

void set_seed_with_primes(Random *p_rnd);

int main(int argc, char *argv[]){

	Random rnd;
        set_seed_with_primes(&rnd);

        double lines_distance = 2;
        double needle_lenght = 1.2;
        double needles_middle_point = 0;
        double needles_inclination = 0;
        double needles_far_left = 0;
        double needles_far_right = 0;

        int number_of_throws = 10E6;
        int number_of_blocks = 1E2;
        int number_of_throws_in_each_block = int(number_of_throws/number_of_blocks);

        double *pi_estimate = new double[number_of_blocks]();
        double *pi2_estimate = new double[number_of_blocks]();
        for(int iblock=0; iblock<number_of_blocks; iblock++){

                double n_hits = 0;
                for(int ithrow=0; ithrow<number_of_throws_in_each_block; ithrow++){

                        needles_middle_point = rnd.Rannyu()*lines_distance;
                        needles_inclination = rnd.Rannyu()*2*M_PI;

                        needles_far_left = needles_middle_point - (needle_lenght/2)*cos(needles_inclination);
                        needles_far_right = needles_middle_point + (needle_lenght/2)*cos(needles_inclination);

                        if( needles_far_left < 0
                            || needles_far_left > lines_distance
                            || needles_far_right < 0
                            || needles_far_right > lines_distance)

                                        n_hits ++;
                }
                pi_estimate[iblock] = (2*needle_lenght*number_of_throws_in_each_block)/(n_hits*lines_distance);
                pi2_estimate[iblock] = pi_estimate[iblock]*pi_estimate[iblock];
        }

        prog_average_std_dev_block_method("results/EX013.dat", pi_estimate, pi2_estimate, number_of_blocks);

	return 0;
}

// Functions definitions

void set_seed_with_primes(Random *p_rnd){
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
				p_rnd->SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
};
