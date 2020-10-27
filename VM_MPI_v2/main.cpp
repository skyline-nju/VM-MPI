#include "config.h"
#include "mpi.h"
#include "communicator2D.h"
#include "domain2D.h"
#include "cellList2D.h"
#include "rand.h"
#include "run2D.h"
#include "particle2D.h"

using namespace std;

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  double Lx = atof(argv[1]);
  double Ly = atof(argv[2]);
  double eta = atof(argv[3]);
  double rho0 = atof(argv[4]);
  double v0 = atof(argv[5]);
  int n_step = atof(argv[6]);
  unsigned long long seed = atoi(argv[7]);
  std::string ini_mode = argv[8];

  Vec_2<double> gl_l(Lx, Ly);
  double eps = 0;
  int seed2 = 0;
  int snap_dt = 10000;
  int field_dt = 10000;
  int field_dx = 4;

  int gl_par_num = int(gl_l.x * gl_l.y * rho0);

  run(gl_par_num, gl_l, eta, eps, v0, n_step, ini_mode,
      seed, seed2, snap_dt, field_dt, field_dx, MPI_COMM_WORLD);
  MPI_Finalize();
}
