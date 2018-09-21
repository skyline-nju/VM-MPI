#include "config.h"
#include "mpi.h"
#include "communicator2D.h"
#include "domain2D.h"
#include "cellList2D.h"
#include "rand.h"
#include "run2D.h"
#include "particle2D.h"
#ifdef OUTPUT_ON
#include "exporter2D.h"
#endif

using namespace std;

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  double L = atoi(argv[1]);
  double eta = atof(argv[2]);
  double h0 = atof(argv[3]);
  int t_half = atoi(argv[4]);
  Vec_2<double> gl_l(L, L);
  double rho0 = 1;
  unsigned long long seed = 1;
  int n_period = 2000;
  double v0 = 0.5;
  double r_cut = 1.;
  int gl_par_num = int(gl_l.x *gl_l.y * rho0);

  run_osc_field(gl_par_num, gl_l, eta, h0, t_half, n_period, seed, r_cut, v0);
  MPI_Finalize();
  return 0;
}