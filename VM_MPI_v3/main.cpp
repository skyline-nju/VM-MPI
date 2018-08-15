#include "mpi.h"
#include "communicator3D.h"
#include "domain3D.h"
#include "cellList3D.h"
#include "particle3D.h"
#include "rand.h"
#include "run.h"
using namespace std;

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  double L = 32;
  Vec_3<double> gl_l(L, L, L);
  double eta = 0.2;
  double rho0 = 1.;
  unsigned long long seed = 12;
  Ranq1 myran(seed);
  int n_step = 1000;
  double v0 = 0.5;
  double r_cut = 1.;
  int gl_par_num = int(gl_l.x * gl_l.y * gl_l.z * rho0);

  typedef UniNode<VicsekPar_3> node_t;
  CellListNode_3<node_t> *cl;
  std::vector<node_t> p_arr;
  Domain_3 dm(gl_l, &cl, r_cut);
  ini_rand(p_arr, gl_par_num, myran, *cl, dm);
  auto interact = [&dm, cl] () {
    cal_force(*cl, dm);
  };
  auto integ = [&dm, cl, &myran, eta, v0](std::vector<node_t> &par_arr) {
   integrate(par_arr, myran, *cl, dm, eta, v0);
  };

  run(p_arr, interact, integ, n_step, 200);
  MPI_Finalize();
}