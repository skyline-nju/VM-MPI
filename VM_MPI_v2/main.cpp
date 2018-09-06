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

#ifdef OUTPUT_ON
  ini_output(gl_par_num, eta, eps, n_step, seed, gl_l, dm.domain_sizes());
  LogExporter log_ex(10000);
  OrderParaExporter phi_ex(100);
  const int field_record_dt = 1000;
  const int cg_box_size = 2;
  FieldExporter field_ex(field_record_dt, 0, cg_box_size,
                         dm.domain_rank(), dm.gl_cells_size(), dm.cells_size());
  auto out = [&log_ex, &phi_ex, &field_ex, my_rank](int i, std::vector<node_t> &par_arr) {
    if (my_rank == 0) {
      log_ex.record(i);
    }
    phi_ex.dump(i, par_arr);
    field_ex.dump(i, par_arr);
  };
#endif

  if (eps == 0) {
    auto integ = [&dm, cl, &myran, eta, v0](std::vector<node_t> &par_arr) {
      integrate(par_arr, myran, *cl, dm, eta, v0);
    };
#ifdef OUTPUT_ON
    run(p_arr, gl_par_num, interact, integ, out, n_step);
#else
    run(p_arr, gl_par_num, interact, integ, n_step, 100);
#endif
  } else {
    RandTorqueArr_2 torques(eps, myran, dm.origin(), dm.cells_size(), dm.gl_cells_size());
    auto integ = [&dm, cl, &myran, eta, v0, &torques](std::vector<node_t> &par_arr) {
      integrate(par_arr, myran, *cl, dm, eta, v0, torques);
    };
#ifdef OUTPUT_ON
    run(p_arr, gl_par_num, interact, integ, out, n_step);
#else
    run(p_arr, gl_par_num, interact, integ, n_step, 100);
#endif
  }

  MPI_Finalize();
  return 0;
}