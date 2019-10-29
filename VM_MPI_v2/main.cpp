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

  double Lx = atof(argv[1]);
  double Ly = atof(argv[2]);
  double eta = atof(argv[3]);
  double eps = atof(argv[4]);
  int seed = atoi(argv[5]);
  int n_step = atoi(argv[6]);
  Vec_2<double> gl_l(Lx, Ly);

  double rho0 = 1;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  //Ranq1 myran(seed + my_rank);
  double v0 = 0.5;
  double r_cut = 1.;
  int gl_par_num = int(gl_l.x *gl_l.y * rho0);

  run_mult_bands(gl_par_num, gl_l, eta, eps, seed, n_step);
}
