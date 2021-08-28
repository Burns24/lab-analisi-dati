#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random_walk.h"
#include "random.h"

using namespace std;

random_walk :: random_walk(){}

random_walk :: ~random_walk(){}

void random_walk :: set_step_lenght(int my_step_lenght){
    /*
    Metodo per impostare la lunghezza del RW.
    */
    step_lenght = my_step_lenght;
};

void random_walk :: set_steps_number(int my_steps_number){
    /*
    Metodo per impostare il numero di step del RW.
    */
    steps_number = my_steps_number;
};

void random_walk :: set_prob_backw(double my_prob_backw){
    /*
    Metodo per impostare la probabilità di uno step all'indietro .
    */
    prob_backw = my_prob_backw;
};

double random_walk :: euclidean_distance(double* pointA, double* pointB){
    /*
    Metodo che restituisce la distanza euclidea date le coordinate di partenza e arrivo.
    */
    return sqrt( pow(pointA[0]-pointB[0],2.) + pow(pointA[1]-pointB[1],2.) + pow(pointA[2]-pointB[2],2.) );
};

void random_walk :: square_lattice(double* random_vec, double* origin, double* walker_head, int repetition, double* coin_toss, double *sum_RW_distance_vec, double *sum_sqr_RW_distance_vec){
    /*
    Metodo che aggiorna i due vettori dati con la somma su tutte le ripetioni della distanza
    percorsa da RW e con la somma su tutte le ripetizioni della distanza al quadrato percorsa
    dal RW su una griglia cubica.
    Il metodo necessita di:
    - double *random_vec: puntatore ad un vettore di double contenente $steps_number numeri casuali distribuiti uniformemente
                          che saranno usati per decidere in che direzione fare lo step;
    - double *origin: puntatore ad un vettore di double con le coordinate dell'origine;
    - double *walker_head: puntatore a vettore contenente le coordinate della testa del cammino;
    - int repetition: intero contenente il numero di volte che il cammino sarà ripetuto;
    - double *coin_toss:  puntatore ad un vettore di double contenente $steps_number numeri casuali distribuiti uniformemente
                          che sarà usato per decidere se lo step vada fatto indietro o avanti;
    - double *sum_RW_distance_vec: puntatore ad un vettore di double lungo $steps_number
                                   contenente le somme delle distanze lungo le traiettorie;
    - double *sum_sqr_RW_distance_vec: puntatore ad un vettore di double lungo $steps_number
                                       contenente le somme al quadrato delle distanze lungo le traiettorie;
    */
    for(int j=0; j < steps_number; j++){                  // Ciclo sugli steps di una singola traiettoria
        int k = j + repetition*steps_number;
        double RW_distance = 0;

        if(random_vec[k]> 0 && random_vec[k]<= 1/3.){               // Seleziono la direzione x
            if(coin_toss[k] <= prob_backw){ walker_head[0] += -step_lenght; }   // Stabilisco se lo step sia indietro o meno
            else { walker_head[0] += step_lenght; }
        } else if(random_vec[k]> 1/3. && random_vec[k]<= 2/3.){     // Seleziono la direzione y
            if(coin_toss[k] <= prob_backw){ walker_head[1] += -step_lenght; }   // Stabilisco se lo step sia indietro o meno
            else { walker_head[1] += step_lenght; }
        } else if(random_vec[k]> 2/3. && random_vec[k]<= 1){        // Seleziono la direzione z
            if(coin_toss[k] <= prob_backw){ walker_head[2] += -step_lenght; }   // Stabilisco se lo step sia indietro o meno
            else { walker_head[2] += step_lenght; }
        }

        RW_distance = euclidean_distance(origin, walker_head);
        sum_RW_distance_vec[j] += RW_distance;
        sum_sqr_RW_distance_vec[j] += pow(RW_distance,2.);
    }
};

void random_walk :: continuum(double* random_theta_vec, double* random_phi_vec, double* origin, double* walker_head, int repetition, double* coin_toss, double *sum_RW_distance_vec, double *sum_sqr_RW_distance_vec){
    /*
    Metodo che aggiorna i due vettori dati con la somma su tutte le ripetioni della distanza
    percorsa da RW e con la somma su tutte le ripetizioni della distanza al quadrato percorsa
    dal RW nel continuo.
    Il metodo necessita di:
    - double* random_theta_vec: puntatore ad un vettore di double contenente $steps_number numeri casuali distribuiti uniformemente
                                in [0, pi] usati per decidere in che direzione theta fare lo step;
    - double* random_phi_vec: puntatore ad un vettore di double contenente $steps_number numeri casuali distribuiti uniformemente
                                in [0, 2pi] usati per decidere in che direzione phi fare lo step;
    - double *origin: puntatore ad un vettore di double con le coordinate dell'origine;
    - double *walker_head: puntatore a vettore contenente le coordinate della testa del cammino;
    - int repetition: intero contenente il numero di volte che il cammino sarà ripetuto;
    - double *coin_toss:  puntatore ad un vettore di double contenente $steps_number numeri casuali distribuiti uniformemente
                          che sarà usato per decidere se lo step vada fatto indietro o avanti;
    - double *sum_RW_distance_vec: puntatore ad un vettore di double lungo $steps_number
                                   contenente le somme delle distanze lungo le traiettorie;
    - double *sum_sqr_RW_distance_vec: puntatore ad un vettore di double lungo $steps_number
                                      contenente le somme al quadrato delle distanze lungo le traiettorie;
    */

    for(int j=0; j < steps_number; j++){          // Ciclo sugli step di una traiettoria
        int k = j + repetition*steps_number;
        double RW_distance = 0;

        if(coin_toss[k] <= prob_backw){           // Decido se fare lo step indietro o avanti
            walker_head[0] += -step_lenght*sin(random_theta_vec[k])*cos(random_phi_vec[k]); // Aggiorno posizione x
            walker_head[1] += -step_lenght*sin(random_theta_vec[k])*sin(random_phi_vec[k]); // Aggiorno posizione y
            walker_head[2] += -step_lenght*cos(random_theta_vec[k]);                        // Aggiorno posizione z
        } else {
            walker_head[0] += step_lenght*sin(random_theta_vec[k])*cos(random_phi_vec[k]);
            walker_head[1] += step_lenght*sin(random_theta_vec[k])*sin(random_phi_vec[k]);
            walker_head[2] += step_lenght*cos(random_theta_vec[k]);
        }

        RW_distance = euclidean_distance(origin, walker_head);
        sum_RW_distance_vec[j] += RW_distance;
        sum_sqr_RW_distance_vec[j] += pow(RW_distance,2.);
    }
};
