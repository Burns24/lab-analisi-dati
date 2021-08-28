/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

/*
If you want to use this code to compute the various observables in function
of the temperature in a certain range, use the bash script `observables_vs_temperature.sh`
and make sure to modify the first row of the `input.dat` file with the string `TEMPERATURE`
and to specify the temperatures to simulate in the file `temperatures.dat`.
*/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  Input(); //Inizialization
  Equilibration(); // Equilibration
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> neqstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Number of equilibration steps = " << neqstep << endl << endl;

  ReadInput >> restart; // if=1 will restart from a configuration
  ReadInput.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

//initial configuration
  if(restart){
      cout << "The simulation will restart from the configuration saved in 'config.final'." << endl;
      ifstream config("config.final");
      if(config.is_open()==0){
          cout << "Error opening the file!" << endl;
      }
      for (int i=0; i<nspin; ++i){
          config >> s[i];
      }
      config.close();
  } else{
      for (int i=0; i<nspin; ++i)
      {
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
      }
  }
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Equilibration(void){
    cout << "Equilibration..." << endl;
    for (int ieqstep=0; ieqstep < neqstep; ieqstep++){
        Move(metro);
    }
}
void Move(int metro)
{
  int o;
  double delta_energy, flip_probability;
  double probability_up;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
        delta_energy = Boltzmann(s[o],o);
        if(delta_energy<=0){
            s[o] = -s[o];
        }
        else{
            flip_probability = exp(-1*beta*delta_energy);
            double r = rnd.Rannyu();
            if (r <= flip_probability){
                s[o] = -s[o];
                accepted ++;
            }
            attempted++;
        }
    }
    else //Gibbs sampling
    {
        probability_up = 1/(1 + exp(-2 * beta *( J * (s[Pbc(o-1)] + s[Pbc(o+1)] ) + h ) ) );
        double r = rnd.Rannyu();
        if (r <= probability_up) s[o] = +1;
        else s[o] = -1;
        accepted ++;
        attempted++;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = 2 * J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) + 2 * h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0;
  double m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
}

void Reset(int iblk) //Reset block averages
{

   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{

   ofstream Ene, Heat, Mag, Chi;
   ofstream Ene_T, Heat_T, Mag_T, Chi_T;
   const int wd=12;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    Ene.open("output.ene.0",ios::app); // Energy
    stima_u = blk_av[iu]/blk_norm/(double)nspin;
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << temp << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Heat.open("output.heat_capacity.0",ios::app); // Heat capacity
    stima_u2 = (blk_av[ic]/blk_norm);
    stima_c = beta*beta*(stima_u2 - (blk_av[iu]/blk_norm)*(blk_av[iu]/blk_norm))/(double)nspin;
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << temp << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

    Mag.open("output.mag.0",ios::app); // Magnetization
    stima_m = blk_av[im]/blk_norm/(double)nspin;
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << temp << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

    Chi.open("output.chi.0",ios::app); // Magnetic susceptibility
    stima_x = beta*(blk_av[ix]/blk_norm)/(double)nspin;
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << temp << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();

    if(iblk==nblk){
        Ene_T.open("results/temp/energy_vs_temperature.dat",ios::app);
        Ene_T << temp << " " << glob_av[iu]/(double)iblk << " " << err_u << endl;
        Ene_T.close();

        Heat_T.open("results/temp/heat_capacity_vs_temperature.dat",ios::app);
        Heat_T << temp << " " << glob_av[ic]/(double)iblk << " " << err_c << endl;
        Heat_T.close();

        Mag_T.open("results/temp/magnetization_vs_temperature.dat",ios::app);
        Mag_T << temp << " " << glob_av[im]/(double)iblk << " " << err_m << endl;
        Mag_T.close();

        Chi_T.open("results/temp/magnetic_susceptibility_vs_temperature.dat",ios::app);
        Chi_T << temp << " " << glob_av[ix]/(double)iblk << " " << err_x << endl;
        Chi_T.close();
    }

// INCLUDE YOUR CODE HERE

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
