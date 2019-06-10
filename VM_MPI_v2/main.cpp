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
  std::cout << "rank =" << my_rank << std::endl;
  //double L = atof(argv[1]);
  double L = 100;
  Vec_2<double> gl_l(L, L);
  //double eta = atof(argv[2]);
  double eta = 0.1;
  double rho0 = 1;
  unsigned long long seed = 1;
  Ranq1 myran(seed + my_rank);
  //int n_step = atoi(argv[4]);
  int n_step = 10000;
  double v0 = 0.5;
  double r_cut = 1.;
  int gl_par_num = int(gl_l.x *gl_l.y * rho0);
  int snap_inteval = 5000;
  int snap_block_size = 2;
  //int box_len_birth = 4;

  run_test(gl_par_num, gl_l, eta, seed, n_step);

}
