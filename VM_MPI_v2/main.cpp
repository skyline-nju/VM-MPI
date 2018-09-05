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
  int my_rank, tot_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  double L = 32;
  Vec_2<double> gl_l(L, L);
  double eta = 0.2;
  double eps = 0;
  double rho0 = 1;
  unsigned long long seed = 1;
  Ranq1 myran(seed + my_rank);
  int n_step = 1000000;
  double v0 = 0.5;
  double r_cut = 1.;
  int gl_par_num = int(gl_l.x *gl_l.y * rho0);

  typedef BiNode<VicsekPar_2> node_t;
  CellListNode_2<node_t> *cl;
  std::vector<node_t> p_arr;
  Domain_2 dm(gl_l, &cl, r_cut);
  ini_rand(p_arr, gl_par_num, myran, *cl, dm);

  // pair-wise interaction
  auto interact = [&dm, cl](std::vector<node_t> &par_arr) {
    cal_force(par_arr, *cl, dm);
  };

  if (eps == 0) {
    auto integ = [&dm, cl, &myran, eta, v0](std::vector<node_t> &par_arr) {
      integrate(par_arr, myran, *cl, dm, eta, v0);
    };
    run(p_arr, gl_par_num, interact, integ, n_step, 100);
  } else {
    RandTorqueArr_2 torques(eps, myran, dm.origin(), dm.cells_size(), dm.gl_cells_size());
    auto integ = [&dm, cl, &myran, eta, v0, &torques](std::vector<node_t> &par_arr) {
      integrate(par_arr, myran, *cl, dm, eta, v0, torques);
    };
    run(p_arr, gl_par_num, interact, integ, n_step, 100);
  }

  MPI_Finalize();
  return 0;
}