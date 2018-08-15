#include "domain2d.h"
#include <iostream>
using namespace std;

void divide_domain(int tot, int &nrows, int &ncols) {
  switch (tot) {
  case 2:
    nrows = 2;
    ncols = 1;
    break;
  case 4:
    nrows = 2;
    ncols = 2;
    break;
  case 6:
    nrows = 3;
    ncols = 2;
    break;
  case 8:
    nrows = 4;
    ncols = 2;
    break;
  case 12:
    nrows = 4;
    ncols = 3;
    break;
  default:
    cout << "np = " << tot << " is wrong!\n";
    exit(1);
  }
}

Domain2d::Domain2d(const cmdline::parser & cmd) {
  eta = cmd.get<double>("eta");
  eps = cmd.get<double>("eps");
  rho0 = cmd.get<double>("rho0");
  gl_Lx = cmd.get<double>("Lx");
  gl_Ly = cmd.get<double>("Ly");
  tot_steps = cmd.get<int>("nstep");
  unsigned long long seed = cmd.get<unsigned long long>("seed");

  MPI_Comm_size(MPI_COMM_WORLD, &tot_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &myrank);
  int nrows_domain;
  int ncols_domain;
  divide_domain(tot_rank, nrows_domain, ncols_domain);
  int my_row_domain = myrank / ncols_domain;
  int my_col_domain = myrank % ncols_domain;

  Lx = gl_Lx / ncols_domain;
  Ly = gl_Ly / nrows_domain;
  x0 = Lx * my_col_domain;
  y0 = Ly * my_row_domain;
  ncols_cell = int(Lx) + 2;
  nrows_cell = int(Ly) + 2;
  cell = NULL;

  particle = NULL;
  nPar = MAX_PAR = end_pos = max_num_comm = 0;

  myran = new Ran(seed + myrank);
  disorder_free = eps > 0 ? false : true;

  // Initialize output
}

Domain2d::~Domain2d() {
  delete[] particle;
  delete[] cell;
}

void Domain2d::create_particle_random(const cmdline::parser & cmd,
                                      int tot_nPar, double magnification) {
  nPar = tot_nPar / tot_rank;
  MAX_PAR = int(nPar * magnification);
  particle = new Node[MAX_PAR];
  end_pos = nPar;
  Node::ini_random(particle, nPar, myran, Lx, Ly, x0, y0);
  max_num_comm = nPar / nrows_cell * 10 * 4;
  
}
