#include "mpi.h"
#include "communicator3D.h"
#include "singleDomain3D.h"
#include "domain3D.h"
#include "cellList3D.h"
using namespace std;

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  double L = 49;
  Vec_3<double> gl_l(L, L, L);
  double eta = 0.2;
  double rho0 = 1.;
  unsigned long long seed = 12;
  Ranq1 myran(seed);
  int n_step = 1000;
  double v0 = 0.5;
  double r_cut = 1.;
  int gl_par_num = int(gl_l.x * gl_l.y * gl_l.z * rho0);

  int neighbor_proc[3][2];
  Vec_3<int> domain_rank;
  Vec_3<bool> flag_ext;

  const Vec_3<int> domain_size = divide_cubic_domain();
  find_neighbor_proc(domain_size, domain_rank, neighbor_proc, flag_ext);

  Vec_3<int> gl_cell_size;
  Vec_3<int> my_cell_size;
  Vec_3<int> my_first_cell;
  Vec_3<double> cell_len;

  divide_cell(gl_l, r_cut, gl_cell_size, cell_len);
  divide_domain(domain_size, domain_rank, gl_cell_size,
                my_cell_size, my_first_cell);

  Vec_3<double> my_l = cell_len * my_cell_size;
  Vec_3<double> my_origin = cell_len * my_first_cell;

  CellListNode_3<UniNode<VicsekPar_3>> cl(my_cell_size, cell_len, my_origin, flag_ext);

  cl.show_origin();
  //std::vector<UniNode<VicsekPar_3>> p_arr;
  //Domain_3 dm(my_l);
  //dm.ini_rand(p_arr, n_par_gl, myran, cl);
  //dm.eval_elapsed_time(p_arr, myran, cl, eta, v0, n_step, 200);

  MPI_Finalize();
};