/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "LSN_Exercise_101.h"

using namespace std;

int main(){
        // Inizialization
        input();
        // Simulation
        cout << "Start the simulation...\n\n";
        while(temperature >= delta_temperature){
                cout << "### Temperature = " << temperature << " ###\n";
                accepted_mutations=0;
                attempted_mutations=0;
                move();
                temperature -= delta_temperature;
        }
        print_path();

        return 0;
}

// Functions
void input(){
        int seed[4];
        int p1, p2;
        ifstream Primes("Primes");
        if (Primes.is_open()){
            Primes >> p1 >> p2 ;
        } else cerr << "PROBLEM: Unable to open Primes" << endl;
        Primes.close();

        ifstream input("seed.in");
        if (input.is_open()){
                    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                    rnd.SetRandom(seed,p1,p2);
                    input.close();
        } else cerr << "PROBLEM: Unable to open seed.in" << endl;

        cout << "Simulated Annealing algorithm to solve\n"
             << "the Traveling Salesman Problem (TSP).   " << endl << endl;

        ifstream ReadInput;
        ReadInput.open("input_exercise_101.dat");
        if (ReadInput.is_open()){
                ReadInput >> initial_temperature;
                ReadInput >> delta_temperature;
                ReadInput >> mc_steps;
        } else cerr << "PROBLEM: Unable to open input_exercise_101.dat" << endl;

        cout << "Parameters used:\n";
        cout << "Initial temperature = " << initial_temperature << endl;
        cout << "Delta temperature = " << delta_temperature << endl;
        cout << "Monte Carlo steps per temperature = " << mc_steps << endl << endl;

        // Read city positions from world_config.dat
        ifstream ReadConfig;
        ReadConfig.open("world_config.dat");
        if (ReadConfig.is_open()){
                ReadConfig >> number_of_cities;
                world.resize(number_of_cities);
                for(int i=0; i<number_of_cities; i++){
                        ReadConfig >> world[i].x;
                        ReadConfig >> world[i].y;
                }
        } else cerr << "PROBLEM: Unable to open world_config.dat" << endl;

        ReadInput.close();
        ReadConfig.close();

        // Initialize variables and path to start the simulation
        iprint = mc_steps/10;
        temperature = initial_temperature;
        path.InitializeChromosome(number_of_cities,&world,&rnd);
};
void move(){
        for(int istep=0; istep < mc_steps; istep++){
                old_path_distance = path.GetPathL1Distance();
                mutate_path();
                new_path_distance = new_path.GetPathL1Distance();

                double probability = exp( (1./temperature)*
                                        (old_path_distance - new_path_distance) );
                if(rnd.Rannyu()<= probability){
                        path = new_path;
                        accepted_mutations ++;
                }
                attempted_mutations ++;

                double path_distance = path.GetPathL1Distance();

                print_distance(path_distance);
                if(istep%iprint==0){
                    cout<<"MC step= "<< istep
                        <<"   Acc. rate= "<< (double)accepted_mutations/attempted_mutations
                        <<"   Path distance= " << path_distance <<endl;
                }
        }
};
void mutate_path(){
        new_path = path;
        int selected_mutation = rnd.Rannyu(0,5);

        if(selected_mutation == 0)
              new_path.PermutePath();
        if(selected_mutation == 1)
              new_path.BlockPermutePath();
        if(selected_mutation == 2)
              new_path.ShiftPath();
        if(selected_mutation == 3)
              new_path.PartialShiftPath();
        if(selected_mutation == 4)
              new_path.InvertPath();
}
void print_distance(double my_path_distance){
    ofstream distance;
    distance.open("results/exercise_10.1/path_distance.dat",ios::app);

    distance << cont << " " << temperature << " " << my_path_distance << endl;

    cont ++;
    distance.close();
}
void print_path(){
    ofstream my_path;
    my_path.open("results/exercise_10.1/path.dat");

    for(int i=0; i<number_of_cities; i++){
            int index=path.GetPath()[i];
            my_path << world[index].x << "   " << world[index].y <<endl;
    }
    my_path << world[path.GetPath()[0]].x << "   " << world[path.GetPath()[0]].y <<endl;

    my_path.close();
}
