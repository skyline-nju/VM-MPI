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
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  double L = atof(argv[1]);
  Vec_2<double> gl_l(L, L);
  double eta = atof(argv[2]);
  double eps = atof(argv[3]);
  double rho0 = 1;
  unsigned long long seed = 1;
  Ranq1 myran(seed + my_rank);
  int n_step = atoi(argv[4]);
  double v0 = 0.5;
  double r_cut = 1.;
  int gl_par_num = int(gl_l.x *gl_l.y * rho0);
  int snap_inteval = 5000;
  int snap_block_size = 2;
  int box_len_birth = 4;

  //run(gl_par_num, gl_l, eta, eps, seed, n_step, n_save_snap, block_size);

  run_birth_death(gl_par_num, gl_l, eta, eps * 1e-8, seed, n_step,
                  snap_inteval, snap_block_size, box_len_birth);
}
