/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "chromosome.h"

using namespace std;

// Constructor
Chromosome :: Chromosome(){};
// Destructor
Chromosome :: ~Chromosome(){};
// Copy constructor
Chromosome :: Chromosome(const Chromosome& my_chromosome){
        this->number_of_cities = my_chromosome.number_of_cities;
        this->p_world = my_chromosome.p_world;
        this->path = my_chromosome.path;
        this->p_rnd = my_chromosome.p_rnd;
};
// Methods
void Chromosome :: InitializeChromosome(int my_number_of_cities, std::vector<city> *my_p_world,
                                        Random *my_p_rnd){
        this->number_of_cities = my_number_of_cities;
        this->p_world = my_p_world;
        this->p_rnd = my_p_rnd;
        // Initialize path with a random permutation
        int number_of_permutation = this->number_of_cities*10;
        path.resize(this->number_of_cities);
        for(int i=0; i < this->number_of_cities; i++)
                this->path[i] = i;
        for(int i=0; i < number_of_permutation; i++)
                this->PermutePath();
        // Check if the created path fulfill the bonds
        this->CheckPath();
};
void Chromosome :: CheckPath(){
        if(this->path[0]!=0){
                cerr<<"Error: path not fulfilling bonds (not starting from the 1^ city)!"<<endl;
                exit (EXIT_FAILURE);
        }
        for(int i=0; i<this->number_of_cities; i++){
                if(this->path[i]<0 || this->path[i]>this->number_of_cities){
                        cerr<<"Error: path not fulfilling bonds (wrong value)!"<<endl;
                        exit (EXIT_FAILURE);
                }
        }
        for(int i=0; i<this->number_of_cities; i++){
                for(int j=i+1; j<this->number_of_cities; j++){
                        if(this->path[i]==this->path[j]){
                                cerr<<"Error: path not fulfilling bonds (duplicated value)!"<<endl;
                                exit (EXIT_FAILURE);
                        }
                }
        }
};
std::vector<int> Chromosome :: GetPath(){
        return this->path;
};
void Chromosome :: SetPath(std::vector<int> my_path){
        this->path = my_path;
};
double Chromosome :: GetPathL1Distance(){
        double l1_distance=0;
        int j=0,k=0;
        for(int i=0; i<(this->number_of_cities-1); i++){
                j = this->path[i];
                k = this->path[i+1];
                l1_distance += sqrt(get_R2_square_norm((*this->p_world)[j].x,(*this->p_world)[k].x,
                                                 (*this->p_world)[j].y, (*this->p_world)[k].y));
        }
        // Add last trip
        j=this->path[number_of_cities-1];
        k=this->path[0];
        l1_distance += sqrt(get_R2_square_norm((*this->p_world)[j].x,(*this->p_world)[k].x,
                                         (*this->p_world)[j].y, (*this->p_world)[k].y));
        return l1_distance;
};
double Chromosome :: GetPathL2Distance(){
        double l2_distance=0;
        int j=0,k=0;
        for(int i=0; i<(this->number_of_cities-1); i++){
                j = this->path[i];
                k = this->path[i+1];
                l2_distance += get_R2_square_norm((*this->p_world)[j].x,(*this->p_world)[k].x,
                                            (*this->p_world)[j].y, (*this->p_world)[k].y);
        }
        // Add last trip
        j=this->path[number_of_cities-1];
        k=this->path[0];
        l2_distance += get_R2_square_norm((*this->p_world)[j].x,(*this->p_world)[k].x,
                                    (*this->p_world)[j].y, (*this->p_world)[k].y);
        return l2_distance;
};
// Mutations
void Chromosome :: PermutePath(){
        // Permute path leaving the first city unchanged
        vector<int> permuted_path=this->path;

        int j=this->p_rnd->Rannyu(1, this->path.size());
        int k=this->p_rnd->Rannyu(1, this->path.size());

        permuted_path[j]=this->path[k];
        permuted_path[k]=this->path[j];

        this->path = permuted_path;
};
void Chromosome :: BlockPermutePath(){
        vector<int> permuted_path=this->path;
        int half_size=this->path.size()/2.;
        int begin=this->p_rnd->Rannyu(1, half_size);
        int end=this->p_rnd->Rannyu(begin, half_size);

        for(int i=begin; i<end; i++){
                permuted_path[i]=this->path[i+half_size];
                permuted_path[i+half_size]=this->path[i];
        }

        this->path = permuted_path;
};
void Chromosome :: ShiftPath(){
        vector<int> shifted_path=this->path;
        int index=p_rnd->Rannyu(1, this->path.size());

        rotate(shifted_path.begin()+1, shifted_path.begin()+index, shifted_path.end());

        this->path = shifted_path;
};
void Chromosome :: PartialShiftPath(){
        vector<int> shifted_path=this->path;
        int begin=this->p_rnd->Rannyu(1, this->path.size()/2.);
        int end=this->p_rnd->Rannyu(begin, this->path.size());
        int index=this->p_rnd->Rannyu(begin, end);

        rotate(shifted_path.begin()+begin, shifted_path.begin()+index ,shifted_path.begin()+end);

        this->path = shifted_path;
};
void Chromosome :: InvertPath(){
        vector<int> inverted_path=this->path;
        int begin=this->p_rnd->Rannyu(1, this->path.size());
        int end=this->p_rnd->Rannyu(begin, this->path.size());

        for(int i=0; i<end-begin; i++){
                inverted_path[i+begin]=this->path[end-i-1];
        }

        this->path = inverted_path;
};
// Functions
double get_R2_square_norm(double x1, double x2, double y1, double y2){
        return (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
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
