/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __LSN_EXERCISE_102__
#define __LSN_EXERCISE_102__

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"
#include "city.h"
#include "chromosome.h"

// Variables
Random rnd;
int mc_steps=0;
int cont=0;
int iprint=0;
// World variables
int number_of_cities=0;
std::vector<city> world;
// Simulated Annealing
double temperature=0, initial_temperature=0;
double delta_temperature=0;
int attempted_mutations=0, accepted_mutations=0;
double old_path_distance=0, new_path_distance=0;
Chromosome path, new_path;

// Functions
void input(void);
void move(void);
void mutate_path(void);
void print_distance(double);
void print_path(void);

#endif
