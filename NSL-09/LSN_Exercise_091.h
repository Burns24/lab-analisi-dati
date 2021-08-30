/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __LSN_EXERCISE_091__
#define __LSN_EXERCISE_091__

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"
#include "city.h"
#include "chromosome.h"
#include "population.h"

//### Variables
int nprint=0;
// Random numbers
Random rnd;

// World variables
int number_of_cities;
std::vector<city> world;

// Genetic Algorithm
int number_of_generations=0;
int population_size=0;
Population path_population;
std::vector<Chromosome> new_path_population, offspring(2);
double r=0;
double permutation_probability=0;
double block_permutation_probability=0;
double shift_probability=0;
double partial_shift_probability=0;
double inversion_probability=0;
double crossover_probability=0;

//### Functions
void Input(void);
int RiggedRoulette(double, int);
void PrintPathL1Distances(Population, int);
void PrintBestPath(Population);
std::vector<Chromosome> GetMutatedChromosomes(Chromosome,Chromosome);

#endif
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
