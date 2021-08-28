#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"

using namespace std;

int main (int argc, char *argv[]){

  Random rnd;  //inizializzo variabile Random
  int seed[4];
  int p1, p2;  //paramentri di input per generatore

// setto i parametri da file seed.in e Primes

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

// genero i numeri casuali e setto i parametri per la media a blocchi

  int M = 1E6;
	int N = 1E2;
	int L = static_cast<int>(M/N);
	double *random_vec = new double[M]();

	for(int i=0; i < M; i++){	
		random_vec[i] = rnd.Rannyu();
	}

//stimo integrale I con distribuzione uniforme


	double *average1 = new double[N]();		// definisco il vettore delle medie
	double *average_sqr1 = new double[N]();		// Definisco il vettore delle medie quadratiche
	for(int i=0; i < N; i++){       // Calcolo la media delle mie osservabili e il quadrato delle medie per il calcolo della varianza
		double sum = 0;
		for(int j=0; j < L; j++){
			int k = j + i*L;
			sum +=  M_PI/2*cos(M_PI/2*random_vec[k]);
		}
		average1[i] = sum/L;
		average_sqr1[i] = pow(average1[i],2);
	}
	prog_average_std_dev_block_method("results/EX02_1(1).dat", average1, average_sqr1, N);

  // Stima usando w(x)=2*(1-x) : soluzione semplice con distribuzione facile da generare (con inverso della cumulativa)
          double *average2 = new double[N]();
          double *average_sqr2 = new double[N]();

          for(int i=0; i < N; i++){
                  double sum = 0;
                  for(int j=0; j < L; j++){
                          int k = j + i*L;
                          double importance_x = 1 - sqrt(1-random_vec[k]);
                          sum += M_PI/4*cos(M_PI/2*importance_x)/(1-importance_x);
                  }
                  average2[i] = sum/L;
                  average_sqr2[i] = pow(average2[i],2);
 		 cout << average2[i] << endl;
          }

          prog_average_std_dev_block_method("results/EX02_1(2).dat", average2, average_sqr2, N);


  return 0;
}
