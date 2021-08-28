#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistical_functions.h"
#include "LSN_Exercise_051.h"

using namespace std;

int main(int argc, char *argv[]){

    input();                // Legge il file input.dat e inizializza le variabili come indicato
    set_seed_with_primes(); // Imposta il seed per la generazione di variabili causali prendendo i numeri da primes

//  compute_average_radius_100_state();
  compute_average_radius_210_state();

// Funzione per confrontare cosa succe se si parte da lotano o da vicino, i parametri vanno impostati manualmente nell'implementazione sottostante
//    test_start_far_from_origin();

	return 0;
}

//Implementazione funzioni

void set_seed_with_primes(void){
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

};

void input(){
  ifstream ReadInput;

  cout << "Hydrogen atom           " << endl;
  cout << "Monte Carlo simulation  " << endl << endl;

  ReadInput.open("input.dat");

  ReadInput >> delta;
  ReadInput >> blocks_number;
  ReadInput >> steps_per_block;
  ReadInput >> equilibration_steps;
  ReadInput >> uniform;
  ReadInput >> gauss;

  ReadInput.close();

  if (uniform) cout << "The program perform Metropolis moves with uniform translations." << endl;
  if (gauss) cout << "The program perform Metropolis moves with gaussian translations." << endl;
  cout << "Random step amplitude = " << delta << endl;
  cout << "Number of blocks = " << blocks_number << endl;
  cout << "Number of steps in one block = " << steps_per_block << endl;
  cout << "Equilibration steps = " << equilibration_steps << endl << endl;
  samples_total_number = blocks_number*steps_per_block;

};

// Funzione per la simulazione del raggio dello stato fondamentale
void compute_average_radius_100_state(void){
    cout << "Computing the average radius of the (1,0,0) state" << endl;
    ofstream outfile;
    outfile.open("results/sampled_points_100_state_gauss.dat");
    state_100 = 1;
    for(int iblock=0; iblock < blocks_number; iblock++){
        if(iblock%10==0) cout << "Block number: " << iblock << endl;

        set_random_walk_head(1, 1, 1); // Inizio il random walk nella posizione (1,1,1)
        equilibrate();                 // Equilibro il Random Walk

        acceptance_rate = 0;
        sum_radius = 0;
        for(int ithrow=0; ithrow < steps_per_block; ithrow++){
            move();
            sum_radius += sqrt(x_current*x_current + y_current*y_current + z_current*z_current);
            outfile << x_current << " " << y_current << " " << z_current << endl;
        }

        average_radius[iblock] = sum_radius/steps_per_block;
        average_radius2[iblock] = average_radius[iblock]*average_radius[iblock];
        if(iblock%10==0) cout << "Block acceptance rate: " <<  acceptance_rate/steps_per_block << endl << endl;
    }

    prog_average_std_dev_block_method("results/average_radius_100_state_gauss.dat",
        average_radius, average_radius2, blocks_number);
    outfile.close();
};

// Stato eccitato (2,1,0)
void compute_average_radius_210_state(void){
    cout << "Computing the average radius of the (2,1,0) state" << endl;
    ofstream outfile;
    outfile.open("results/sampled_points_210_state_gauss.dat");
    state_210 = 1;
    for(int iblock=0; iblock < blocks_number; iblock++){
        if(iblock%10==0) cout << "Block number: " << iblock << endl;

        set_random_walk_head(1, 1, 1); // Inizio il random walk nella posizione (1,1,1)
        equilibrate();                 // Equilibrio il random walk

        acceptance_rate = 0;
        sum_radius = 0;
        for(int ithrow=0; ithrow < steps_per_block; ithrow++){
            move();
            sum_radius += sqrt(x_current*x_current + y_current*y_current + z_current*z_current);
            outfile << x_current << " " << y_current << " " << z_current << endl;
        }

        average_radius[iblock] = sum_radius/steps_per_block;
        average_radius2[iblock] = average_radius[iblock]*average_radius[iblock];
        if(iblock%10==0) cout << "Block acceptance rate: " <<  acceptance_rate/steps_per_block << endl << endl;
    }

    prog_average_std_dev_block_method("results/average_radius_210_state_gauss.dat",
        average_radius, average_radius2, blocks_number);
    outfile.close();
};


void test_start_far_from_origin(void){
    cout << "Testing what happens when starting the random walk far from the origin." << endl;
  //Cambio il nome del file di output, la posizione di partenza e lo stato selezionato in base a cosa mi occorre
    ofstream out;
    out.open("results/test_start_far_from_origin_100_state.dat");

    set_random_walk_head(20, 20, 20);
    state_100 = 1;
    for(int istep=0; istep< steps_per_block; istep++){
        move();
        out << sqrt(x_current*x_current + y_current*y_current + z_current*z_current) << endl;
    }
    out.close();
};


// Funzione che setta la posizione attuale del random walk
void set_random_walk_head(double my_x,double my_y,double my_z){
    x_current = my_x;
    y_current = my_y;
    z_current = my_z;
};

// Funzione per equilibrare il random walk
void equilibrate(void){
    for(int ieq=0; ieq< equilibration_steps; ieq++){
        move();
    }
    acceptance_rate = 0;
};


// Funzione per la mossa in 3D
void move(void){

//Calcolo l'ampiezza di probabilità attuale in base allo stato (100 o 210)
    if(state_100) current_probability_amplitude =
        get_probability_amplitude_100_state(x_current, y_current, z_current);
    if(state_210) current_probability_amplitude =
        get_probability_amplitude_210_state(x_current, y_current, z_current);

    if (uniform) make_uniform_step();
    if (gauss) make_gaussian_step();

    if(state_100) new_probability_amplitude =
        get_probability_amplitude_100_state(x_new, y_new, z_new);
    if(state_210) new_probability_amplitude =
        get_probability_amplitude_210_state(x_new, y_new, z_new);

    acceptance = get_acceptance(current_probability_amplitude,
                                new_probability_amplitude);
    double r = rnd.Rannyu();
    if (r <= acceptance){
        set_random_walk_head(x_new, y_new, z_new);
        acceptance_rate ++;
    }
};


double get_probability_amplitude_100_state(double my_x,double my_y,double my_z){
    double probability = (1/(pow(bohr_radius,3)*M_PI))
                         *exp(-2*sqrt(my_x*my_x+my_y*my_y+my_z*my_z)/bohr_radius);
    return probability;
};

double get_probability_amplitude_210_state(double my_x,double my_y,double my_z){
    double probability = (1/(32*pow(bohr_radius,5)*M_PI))*my_z*my_z
                          *exp(-1*sqrt(my_x*my_x+my_y*my_y+my_z*my_z)/bohr_radius);
    return probability;
};

// Step uniforme
void make_uniform_step(void){
    x_new = x_current + delta*(rnd.Rannyu() - 0.5);
    y_new = y_current + delta*(rnd.Rannyu() - 0.5);
    z_new = z_current + delta*(rnd.Rannyu() - 0.5);
};

// Step Gaussiano
void make_gaussian_step(void){
    x_new = x_current + rnd.Gauss(0,delta);
    y_new = y_current + rnd.Gauss(0,delta);
    z_new = z_current + rnd.Gauss(0,delta);
};

// Funzione di accettazione
double get_acceptance(double prob_current, double prob_new){
    double a = prob_new/prob_current;
    if (a < 1) return a;
    else return 1;
};
