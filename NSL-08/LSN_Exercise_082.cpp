/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "LSN_Exercise_082.h"

using namespace std;

int main()
{
  // Inizialization
  Input();
  // Optimize mu and sigma
  if(simulated_annealing==1) OptimizeParameters();
  // Equilibration
  cout << "Equilibration..." << endl << endl;
  for(int ieq=0; ieq<neqstep; ieq++)
      Move(mu, sigma);

  // Simulation
  for(int iblk=1; iblk <= nblk; ++iblk)
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(mu, sigma);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }

  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadParam;

  cout << "Single 1D quantum particle         " << endl;
  cout << "Variational Monte Carlo simulation." << endl << endl;
  cout << "External confinement potential: v(x) = x^4 - 5/2 * x^2" << endl;
  cout << "Approximate ground state function:" << endl;
  cout << "\texp( -(x-mu)**2 / 2*sigma**2 ) + exp( -(x+mu)**2 / 2*sigma**2 )" << endl << endl;


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
  ReadParam.open("parameters.dat");

  ReadInput >> x_pos;
  ReadInput >> delta;
  ReadInput >> neqstep;
  ReadInput >> nblk;
  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Initial x-position = " << x_pos << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of equilibration steps = " << neqstep << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  ReadInput >> simulated_annealing;
  if(simulated_annealing==0)
  {
      ReadParam >> mu;
      ReadParam >> sigma;
      cout << "Simulation will be performed with parameters mu = " << mu
           << " and sigma = " << sigma << endl << endl;
  }
  if(simulated_annealing==1)
  {
      ReadInput >> starting_temp;
      ReadParam >> starting_mu;
      ReadParam >> starting_sigma;

      cout << "In order to choose mu and sigma a Simulated Annealing" << endl;
      cout << "procedure is performed with:                         " << endl;
      cout << "Starting temperature = " << starting_temp << endl;
      cout << "Starting mu = " << starting_mu << endl;
      cout << "Starting sigma = " << starting_sigma << endl << endl;
  }

  ReadInput.close();
  ReadParam.close();

//Prepare arrays for measurements
  ih = 0; // Hamiltonian
  n_props = 1; //Number of observables

}

void OptimizeParameters(void)
{
    /* Function to optimize the parameters of the variational wave function
    (mu and sigma) according to a simulated annealing procedure. Values for which
    <H> is minimized are saved in `opt_parameters.dat`.
    */
    double delta_temp = 0.0005; //funziona: 0.00005
    double delta_opt = 0.01;
    double x0 = 0;
    int noptstep = 1E2, noptcycle=3;
    int nave_ham = 1E4;
    int nprint = 20; // how often to write on file mu and sigma values
    bool opt_verbose = 1; // set to 1 if you want to print temperature, acceptance and other variables during optimization
    int wd = 10;

    cout << "Parameters will be optimized starting every time from x0 = "<< x0 << " and with: " << endl;
    cout << "Delta temperature = " << delta_temp << endl;
    cout << "Cycles of optimization = " << noptcycle << endl;
    cout << "Optimization steps per temperature = " << noptstep << endl;
    cout << "Samples for <H> per optimization step = " << nave_ham << endl;

    ofstream sampled_mu_sigma;
    sampled_mu_sigma.open("results/sampled_mu_sigma.dat", ios::app);
    int opt_accepted=0, opt_attempted=0;
    double h_ave = 0;
    double old_ave_ham=0, new_ave_ham=0;
    double new_mu=0, new_sigma=0;
    mu = starting_mu;
    sigma = starting_sigma;
    for(int ioptcycle=0; ioptcycle<noptcycle; ioptcycle++)
    {
        cout << "\n### Opt. cycle number: " << ioptcycle << endl;

        double temperature = starting_temp;
        while(temperature >= delta_temp)
        {
            opt_accepted = 0;
            opt_attempted = 0;
            for(int ioptstep=0; ioptstep<noptstep; ioptstep++)
            {
                x_pos = x0;
                old_ave_ham = 0;
                for(int i=0; i<nave_ham; i++)
                {
                    Move(mu, sigma);
                    old_ave_ham += Hamiltonian(x_pos, mu, sigma);
                }
                old_ave_ham /= nave_ham;

                new_mu = mu + delta_opt*(rnd.Rannyu() - 0.5);
                new_sigma = sigma + delta_opt*(rnd.Rannyu() - 0.5);

                x_pos = x0;
                new_ave_ham = 0;
                for(int i=0; i<nave_ham; i++)
                {
                    Move(new_mu, new_sigma);
                    new_ave_ham += Hamiltonian(x_pos, new_mu, new_sigma);
                }
                new_ave_ham /= nave_ham;

                double p = exp((1/temperature) * (old_ave_ham-new_ave_ham));
                if(rnd.Rannyu() <= p)
                {
                    mu = new_mu;
                    sigma = new_sigma;
                    opt_accepted ++;
                }
                opt_attempted ++;

                if(ioptstep%nprint==0)
                {
                    x_pos = x0;
                    h_ave=0;
                    for(int i=0; i<nave_ham; i++){
                        Move(mu, sigma);
                        h_ave+=Hamiltonian(x_pos, mu, sigma);
                    }
                    h_ave /= nave_ham;
                    if(opt_verbose==1)
                    {
                        cout <<setw(wd)<< "Temperature = " <<setw(7)<< temperature <<setw(wd) << " Opt. step = "
                             <<setw(wd)<< ioptstep <<setw(wd)<< " Acceptance = "
                             <<setw(3)<< opt_accepted/(double) opt_attempted <<setw(wd)<< " mu = " <<setw(5)<< mu
                             <<setw(wd)<< " sigma = " <<setw(5)<< sigma <<setw(7)<< " <H> = "
                             <<setw(wd)<< h_ave << endl;
                    }
                    sampled_mu_sigma << mu << " " << sigma << " " << h_ave << endl;
                }
                }

            temperature = temperature - delta_temp;
        }
    }
    sampled_mu_sigma.close();

    // Select mu and sigma for which <H> is minimized
    int nsample = noptcycle*(noptstep/nprint)*(starting_temp/delta_temp);
    double *sampled_mu = new double[nsample]();
    double *sampled_sigma = new double[nsample]();
    double *sampled_ham = new double[nsample]();

    ifstream sampled_param;
    sampled_param.open("results/sampled_mu_sigma.dat");
    for(int isample=0; isample<nsample; isample++)
        sampled_param >> sampled_mu[isample] >> sampled_sigma[isample] >> sampled_ham[isample];

    int h_min_index = 0;
    double h_min = sampled_ham[0];
    for(int isample=0; isample<nsample; isample++)
    {
            if(h_min > sampled_ham[isample])
            {
                h_min = sampled_ham[isample];
                h_min_index = isample;
            }
    }

    // Saving optimized parameters
    ofstream opt_parameters;
    opt_parameters.open("parameters.dat");
    opt_parameters << sampled_mu[h_min_index] << "\n" << sampled_sigma[h_min_index];
    opt_parameters.close();

    mu = sampled_mu[h_min_index];
    sigma = sampled_sigma[h_min_index];
    cout << "\nSimulation will be performed with optimized parameters mu = " << mu
         << " and sigma = " << sigma << endl;
}

void Move(double my_mu, double my_sigma)
{
    double old_prob=0, new_prob=0;
    double new_x_pos=0;

    old_prob = VariationalPsi(x_pos,my_mu,my_sigma)*VariationalPsi(x_pos,my_mu,my_sigma);
    new_x_pos = x_pos + delta*(rnd.Rannyu() - 0.5);
    new_prob = VariationalPsi(new_x_pos,my_mu,my_sigma)*VariationalPsi(new_x_pos,my_mu,my_sigma);

    if( rnd.Rannyu() <= new_prob/old_prob )
    {
        x_pos = new_x_pos;
        accepted ++;
    }
    attempted ++;
}

void Measure()
{
    walker[ih] = Hamiltonian(x_pos, mu, sigma);

    ofstream position;
    position.open("results/x_pos.dat",ios::app);
    position << x_pos << endl;
    position.close();

}

double VariationalPsi(double my_x,double my_mu,double my_sigma)
{
    double diff = my_x - my_mu;
    double sum = my_x + my_mu;
    double sigma2 = my_sigma*my_sigma;

    return exp(-1*diff*diff/ (2*sigma2)) + exp(-1*sum*sum / (2*sigma2));
}

double Hamiltonian(double my_x, double my_mu, double my_sigma)
{
    double kin_energy = -0.5 *(
                   exp(-(my_x-my_mu)*(my_x-my_mu)/(2.*my_sigma*my_sigma)) * (my_mu*my_mu-my_sigma*my_sigma+my_x*my_x-2.*my_mu*my_x)/(my_sigma*my_sigma*my_sigma*my_sigma) +
                   exp(-(my_x+my_mu)*(my_x+my_mu)/(2.*my_sigma*my_sigma)) * (my_mu*my_mu-my_sigma*my_sigma+my_x*my_x+2.*my_mu*my_x)/(my_sigma*my_sigma*my_sigma*my_sigma) );

    double pot_energy = my_x*my_x*my_x*my_x - 5./2.*my_x*my_x;
    return kin_energy + pot_energy;
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

    ofstream H_ave;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    H_ave.open("results/h_ave.dat",ios::app);

    stima_ham = blk_av[ih]/blk_norm;
    glob_av[ih] += stima_ham;
    glob_av2[ih] += stima_ham*stima_ham;
    err_ham=Error(glob_av[ih],glob_av2[ih],iblk);

// Hamiltonian
    H_ave << " "<< iblk << " " << glob_av[ih]/(double)iblk << " "<< err_ham << endl;

    cout << "----------------------------" << endl << endl;

    H_ave.close();

}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
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
