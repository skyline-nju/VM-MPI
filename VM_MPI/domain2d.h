#ifndef DOMAIN_2D_H
#define DOMAIN_2D_H
#include <vector>
#include <stack>
#include "mpi.h"
#include "rand.h"
#include "cell.h"
#include "cmdline.h"

class Domain2d
{
public:
  Domain2d(const cmdline::parser &cmd);
  ~Domain2d();

  void create_particle_random(const cmdline::parser &cmd,
                              int tot_nPar, double magnification);
private:
  int tot_rank;
  int myrank;
  int left_rank;
  int right_rank;
  int up_rank;
  int down_rank;

  double eta;
  double eps;
  double rho0;
  double gl_Lx;         // Length of global system in x direction.
  double gl_Ly;
  int tot_steps;

  double x0;
  double y0;
  double Lx;
  double Ly;
  int ncols_cell;
  int nrows_cell;
  Cell *cell;

  Node *particle;
  int nPar;           // Total number of valid particles.
  int MAX_PAR;        // Max number of particles.
  int end_pos;        
  std::stack <Node *> empty_pos;
  int max_num_comm;  // Max number of double float values need to send or receive per communication.
  int nPar_per_proc; // Mean number of particles per process.

  Ran *myran;
  bool disorder_free;
 
  std::time_t beg_time;
  std::ofstream fout; // Output order parameters by proc 0 and log by proc 1.
};
#endif
