#include "mpi.h"
#include "communicator3D.h"
#include "domain3D.h"
#include "cellList3D.h"
#include "particle3D.h"
#include "rand.h"

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

  const auto domain_size = Domain_3::partitioning();
  ExtDomain_3 dm(domain_size, gl_l, r_cut, gl_par_num);
  typedef UniNode<VicsekPar_3_w_flag> node_t;
  std::vector<node_t> p_arr;
  CellListNode_3<node_t> cl(dm.cell_size(), dm.cell_len(), gl_l,
                           dm.origin(), dm.flag_ext());
  //cl.show_para();
  dm.ini_rand(p_arr, myran, cl);
  auto interact = [&dm, &cl] (std::vector<node_t> &par_arr) {
   dm.cal_pair_force(par_arr, cl);
  };
  auto integrate = [&dm, &cl, &myran, eta, v0](std::vector<node_t> &par_arr) {
   dm.integrate(par_arr, myran, cl, eta, v0);
  };

  run(p_arr, interact, integrate, n_step, 200);
  MPI_Finalize();
}