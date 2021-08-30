/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "population.h"

using namespace std;

// Constructor
Population :: Population(){};
// Destructor
Population :: ~Population(){};
// Copy constructor
Population :: Population(const Population& my_population){
        this->population_size = my_population.population_size;
        this->p_rnd = my_population.p_rnd;
        this->population = my_population.population;
        this->crossover_probability = my_population.crossover_probability;
};
// Methods
std::vector<Chromosome> Population :: GetPopulation(){
        return population;
};
void Population :: SetPopulation(std::vector<Chromosome> my_population){
        this->population = my_population;
};
void Population :: InitializePopulation(int my_population_size, int my_number_of_cities,
                                        std::vector<city> *my_p_world, Random *my_p_rnd){
        this->population_size = my_population_size;
        this->p_rnd = my_p_rnd;
        // Initialize population with chromosomes that have
        // uniform distributed path permutation
        this->population.resize(this->population_size);
        for(int i=0; i< this->population_size; i++)
                this->population[i].InitializeChromosome(my_number_of_cities,
                                                   my_p_world,p_rnd);
};
void Population :: SortPopulation(){
        vector<double> distances(population_size);
        vector<int> sorted_indexes(population_size);
        vector<Chromosome> new_population(population_size);
        // Compute L1 distances
        for(int i=0; i<this->population_size; i++){
                distances[i]= this->population[i].GetPathL1Distance();
        }
        // Sort indexes according to the L1 distance of the single path
        sorted_indexes = get_sorted_indexes_vector(distances);
        // Order the population
        for(int i=0; i<this->population_size; i++){
                new_population[i] = this->population[sorted_indexes[i]];
        }
        this->population = new_population;
};
void Population :: SetCrossoverProbability(double my_crossover_probability){
        this-> crossover_probability = my_crossover_probability;
};
std::vector<Chromosome> Population :: GetCrossoveredChromosomes(int k, int l){
        vector<Chromosome> offspring(2);
        vector<vector<int>> offspring_path(2);
        offspring[0] = this->population[k];
        offspring[1] = this->population[l];
        offspring_path[0]=offspring[0].GetPath();
        offspring_path[1]=offspring[1].GetPath();

        if(this->p_rnd->Rannyu()<crossover_probability){
                int index=this->p_rnd->Rannyu(1, offspring_path[0].size());

                vector<int> sequence1=offspring_path[1];
                vector<int> sequence2=offspring_path[0];
                for(unsigned int i=0; i<sequence1.size(); i++){
                        if( find(offspring_path[0].begin()+index, offspring_path[0].end(), sequence1[i])
                                 == offspring_path[0].end() ){
                        sequence1.erase(sequence1.begin()+i);
                        i--;
                        }
                }
                for(unsigned int i=0; i<sequence2.size(); i++){
                        if( find(offspring_path[1].begin()+index, offspring_path[1].end(), sequence2[i])
                                 == offspring_path[1].end() ){
                        sequence2.erase(sequence2.begin()+i);
                        i--;
                        }
                }
                for(unsigned int i=0; i<sequence1.size(); i++){
                        offspring_path[0][i+index]=sequence1[i];
                        offspring_path[1][i+index]=sequence2[i];
                }

                offspring[0].SetPath(offspring_path[0]);
                offspring[1].SetPath(offspring_path[1]);
                // Check if offsprings fullfil the bonds
                offspring[0].CheckPath();
                offspring[1].CheckPath();

               return offspring;
       }else{
               return offspring;
       };

};
// Functions
vector<int> get_sorted_indexes_vector(const vector<double> &v){
        // Initialize original index locations
        vector<int> idx(v.size());
        iota(idx.begin(), idx.end(), 0);
        // Sort indexes based on comparing values in v
        sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

        return idx;
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
