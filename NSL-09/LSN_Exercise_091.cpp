/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "LSN_Exercise_091.h"

using namespace std;

int main()
{
        // Inizialize
        Input();
        // Start simulation
        for(int igeneration=1; igeneration<=number_of_generations; igeneration++){
                if(igeneration%nprint==0)
                        cout<<"Generation number: "<<igeneration<<endl;
                path_population.SortPopulation(); // Order population for selection
                for(int j=0; j<population_size; j=j+2){
                        int k=RiggedRoulette(r, population_size);
                        int l=RiggedRoulette(r, population_size);
                        offspring = path_population.GetCrossoveredChromosomes(k, l);
                        offspring = GetMutatedChromosomes(offspring[0],offspring[1]);
                        new_path_population[j] = offspring[0];
                        new_path_population[j+1] = offspring[1];
                };
                path_population.SetPopulation(new_path_population);
                PrintPathL1Distances(path_population,igeneration);
        };
        PrintBestPath(path_population);

        return 0;
}

//### Functions
void Input(void){
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

        cout << "Solving the Traveling Salesman Problem (TSP)" << endl;
        cout << "via a Genetic Algorithm.                    " << endl << endl;

        ifstream ReadInput, ReadConfig;
        //Read informations from input.dat
        ReadInput.open("input.dat");
        ReadInput >> number_of_generations;
        ReadInput >> population_size;
        ReadInput >> r;
        ReadInput >> permutation_probability;
        ReadInput >> block_permutation_probability;
        ReadInput >> shift_probability;
        ReadInput >> partial_shift_probability;
        ReadInput >> inversion_probability;
        ReadInput >> crossover_probability;
        // Read city positions from world_config.dat
        ReadConfig.open("world_config.dat");
        ReadConfig >> number_of_cities;
        world.resize(number_of_cities);
        for(int i=0; i<number_of_cities; i++){
                ReadConfig >> world[i].x;
                ReadConfig >> world[i].y;
        }

        cout << "Simulation parameters: " << endl;
        cout << "Number of generations = " << number_of_generations << endl;
        cout << "Population size = " << population_size << endl;
        cout << "Number of cities = " << number_of_cities << endl<<endl;

        ReadInput.close();
        ReadConfig.close();

        nprint = number_of_generations/10;

        // Initialize empty population vector
        new_path_population.resize(population_size);
        // Initialize population class
        path_population.InitializePopulation(population_size, number_of_cities,
                                             &world, &rnd);
        path_population.SetCrossoverProbability(crossover_probability);
};
// Selection
int RiggedRoulette(double r, int population_size){
    return pow(rnd.Rannyu(), r)*population_size;
};
void PrintPathL1Distances(Population my_path_population, int igeneration){
        Population sorted_population;
        ofstream best_path, average_best_half_path;
        best_path.open("results/distance_of_best_path.dat",ios::app);
        average_best_half_path.open("results/average_distance_over_best_half_paths.dat",ios::app);

        sorted_population = my_path_population;
        sorted_population.SortPopulation();

        best_path << igeneration << " " << sorted_population.GetPopulation()[0].GetPathL1Distance() <<endl;

        double avg_dist=0;
        for(unsigned int i=0; i<sorted_population.GetPopulation().size()/2; i++){
                avg_dist += sorted_population.GetPopulation()[i].GetPathL1Distance();
        }
        avg_dist /= sorted_population.GetPopulation().size()/2.;

        average_best_half_path << igeneration << " " << avg_dist <<endl;

        best_path.close();
        average_best_half_path.close();
};
void PrintBestPath(Population my_path_population){
        Population sorted_population;
        ofstream best_path;
        best_path.open("results/best_path.dat");

        sorted_population = my_path_population;
        sorted_population.SortPopulation();

        for(int i=0; i<number_of_cities; i++){
                int index=sorted_population.GetPopulation()[0].GetPath()[i];
                best_path << world[index].x << "   " << world[index].y <<endl;
        }
        best_path << world[sorted_population.GetPopulation()[0].GetPath()[0]].x << "   " << world[sorted_population.GetPopulation()[0].GetPath()[0]].y <<endl;

        best_path.close();
};
// Mutation
std::vector<Chromosome> GetMutatedChromosomes(Chromosome my_chromosome1, Chromosome my_chromosome2){
        vector<Chromosome> mutated_chromosome(2);
        mutated_chromosome[0] = my_chromosome1;
        mutated_chromosome[1] = my_chromosome2;

        // Chromosome 1
        if(rnd.Rannyu() < permutation_probability)
              mutated_chromosome[0].PermutePath();
        if(rnd.Rannyu() < block_permutation_probability)
              mutated_chromosome[0].BlockPermutePath();
        if(rnd.Rannyu() < shift_probability)
              mutated_chromosome[0].ShiftPath();
        if(rnd.Rannyu() < partial_shift_probability)
              mutated_chromosome[0].PartialShiftPath();
        if(rnd.Rannyu() < inversion_probability)
                mutated_chromosome[0].InvertPath();
        // Chromosome 2
        if(rnd.Rannyu() < permutation_probability)
              mutated_chromosome[1].PermutePath();
        if(rnd.Rannyu() < block_permutation_probability)
              mutated_chromosome[1].BlockPermutePath();
        if(rnd.Rannyu() < shift_probability)
              mutated_chromosome[1].ShiftPath();
        if(rnd.Rannyu() < partial_shift_probability)
              mutated_chromosome[1].PartialShiftPath();
        if(rnd.Rannyu() < inversion_probability)
                mutated_chromosome[1].InvertPath();

        // Check if mutated paths fullfil the bonds
        mutated_chromosome[0].CheckPath();
        mutated_chromosome[1].CheckPath();

        return mutated_chromosome;
};

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
