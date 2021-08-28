
#ifndef __LSN_Exercise_051__
#define __LSN_Exercise_051__

// generatore numeri casuali
Random rnd;

// Variabili
const double bohr_radius = 1; // In unità di Bohr
double x_new=0, y_new=0, z_new=0; //Posizioni
double x_current=0, y_current=0, z_current=0;
double sum_radius = 0;
double current_probability_amplitude = 0;
double new_probability_amplitude = 0;
double acceptance = 0;
double acceptance_rate = 0;
double delta = 0.5;

// Parametri della simulazione (scelti di default)
int blocks_number = 1E2;
int steps_per_block = 1E4;
int samples_total_number = blocks_number*steps_per_block;
int equilibration_steps = 1E2;
bool state_100=0, state_210=0;
bool uniform =0, gauss=0;

// Blocks statistics
double *average_radius = new double[blocks_number]();
double *average_radius2 = new double[blocks_number]();

// Functions
void input(void);
void set_seed_with_primes(void);
void compute_average_radius_100_state(void);
void compute_average_radius_210_state(void);
void test_start_far_from_origin(void);
void move(void);
void make_uniform_step(void);
void make_gaussian_step(void);
void equilibrate(void);
void set_random_walk_head(double my_x,double my_y,double my_z);
double get_probability_amplitude_100_state(double my_x,double my_y,double my_z);
double get_probability_amplitude_210_state(double x, double y, double z);
double get_acceptance(double p_current, double p_new);

#endif
