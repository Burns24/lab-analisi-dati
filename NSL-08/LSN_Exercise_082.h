/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __LSN_EXERCISE_081__
#define __LSN_EXERCISE_081__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

// simulated annealing
bool simulated_annealing=0;
double starting_temp=0, starting_mu=0, starting_sigma=0;

//parameters, observables
double mu=0.8, sigma=0.6;
double x_pos=0;

const int m_props=10;
int n_props, ih;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_ham, err_ham;

// simulation
int neqstep, nstep, nblk;
double delta;

//functions
void Input(void);
void OptimizeParameters(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(double,double);
void Measure(void);
double VariationalPsi(double,double,double);
double Hamiltonian(double,double,double);
double Error(double,double,int);

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
