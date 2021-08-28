#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"
#include "random_walk.h"

using namespace std;

int main(int argc, char *argv[]){

    // Inizializzazione classe rnd
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


  int N = 1E2;                    // Numero di step
   int L = 1E4;                    // Numero di ripetizioni della traiettoria
   int M = L*N;	        	    // Totale numeri casuali richiesti
   double *coin = new double[M](); // Carico un vettore con numeri casuali distribuiti uniformemente
 for(int i=0; i < M; i++){       // e li userò per decidere se fare lo step avanti o indietro;
       coin[i] = rnd.Rannyu();
   }
   double prob_backw = 0.5;        // Probabilità di fare uno step indietro
   int lattice_constant = 1;       // costante che definisce distanza caratteristica del reticolo
   double *walker_head = new double[3]; // Vettore contenente le coordinate della testa
   double* origin = new double[3];      //Vettore coordinate dell'origine
   origin[0] = 0;
   origin[1] = 0;
   origin[2] = 0;
   random_walk r_walker;
   r_walker.set_step_lenght(lattice_constant);
   r_walker.set_steps_number(N);
   r_walker.set_prob_backw(prob_backw);

   // Stima della distanza 3D percorsa da un Random Walk su di un reticolo cubico
   double *random_vec = new double[M]();   // Definisco e carico un vettore con numeri casuali
 for(int i=0; i < M; i++){	            // distribuiti uniformememnte per decidere la direzione
   random_vec[i] = rnd.Rannyu();       // xyz dello step del random walk;
 }
   double *sum_RW_distance_vec1 = new double[N]();		    // Vettore contente le somme delle distanze su tutte le traiettorie
   double *sum_sqr_RW_distance_vec1 = new double[N]();		// e un vettore delle somme al quadrato

 for(int i=0; i < L; i++){      // Ciclo sulle traiettorie
       walker_head[0] = origin[0];        // Il cammino ha inizio sempre nell'origine
       walker_head[1] = origin[1];
       walker_head[2] = origin[2];
       r_walker.square_lattice(random_vec, origin, walker_head, i, coin, sum_RW_distance_vec1, sum_sqr_RW_distance_vec1);
   }

   ofstream out_file1;                     // Apro un file in cui scrivere i risultati: media e deviazione std.
   out_file1.open("results/EX02_2(1).dat");   // delle distanze percorse dal RW sul reticolo per le varie traiettorie
   for(int j=0; j < N; j++){
      out_file1 << sum_RW_distance_vec1[j]/L << " " << std_dev(sum_RW_distance_vec1[j]/L,sum_sqr_RW_distance_vec1[j]/L,L) << endl;
   }
   out_file1.close();

   // Stima della distanza 3D percorsa da un RW nel continuo
   double theta_min = 0, theta_max = M_PI;
   double *random_theta_vec = new double[M]();		// Definisco e carico un vettore usato per determinare l'angolo theta
 for(int i=0; i < M; i++){	                    // della direzione con numeri distribuiti uniformemente in [0,pi];
   random_theta_vec[i] = rnd.Rannyu()*(theta_max-theta_min) + theta_min;
 }
   double phi_min = 0, phi_max = 2*M_PI;
   double *random_phi_vec = new double[M]();		// Definisco e carico un vettore usato per determinare l'angolo phi
 for(int i=0; i < M; i++){	                    // della direzione con numeri distribuiti uniformemente in [0,2pi];
   random_phi_vec[i] = rnd.Rannyu()*(phi_max-phi_min) + phi_min;
 }
   double *sum_RW_distance_vec2 = new double[N]();		    // Vettore contente le somme delle distanze su tutte le traiettorie
   double *sum_sqr_RW_distance_vec2 = new double[N]();		// e un vettore delle somme al quadrato
   for(int i=0; i < L; i++){              // Ciclo sulle traiettorie
       walker_head[0] = origin[0];        // Il cammino ha inizio sempre nell'origine
       walker_head[1] = origin[1];
       walker_head[2] = origin[2];
       r_walker.continuum(random_theta_vec,random_phi_vec, origin, walker_head, i, coin, sum_RW_distance_vec2, sum_sqr_RW_distance_vec2);
   }
   ofstream out_file2;                     // Apro un file in cui scrivere i risultati: media e deviazione std.
   out_file2.open("results/EX02_2(2).dat");   // delle distanze percorse dal RW sul reticolo per le varie traiettorie
   for(int j=0; j < N; j++){
      out_file2 << sum_RW_distance_vec2[j]/L << " " << std_dev(sum_RW_distance_vec2[j]/L,sum_sqr_RW_distance_vec2[j]/L,L) << endl;
  }
   out_file2.close();

   return 0;

 }
