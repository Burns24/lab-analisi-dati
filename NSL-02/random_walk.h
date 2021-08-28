#ifndef __random_walk__
#define __random_walk__

class random_walk {

private:
    int step_lenght, steps_number;
    double prob_backw;

protected:

public:
  // costruttore
  random_walk();
  // distruttore
  ~random_walk();
  // metodi
  void set_step_lenght(int );
  void set_steps_number(int );
  void set_prob_backw(double );

  double euclidean_distance(double* pointA, double* pointB);

  void square_lattice(double* random_vec, double* origin, double* walker_head, int repetition, double* coin_toss, double *sum_RW_distance_vec, double *sum_sqr_RW_distance_vec);
  void continuum(double* random_theta_vec, double* random_phi_vec, double* origin, double* walker_head, int repetition, double* coin_toss, double *sum_RW_distance_vec, double *sum_sqr_RW_distance_vec);
};

#endif // __random_walk_
