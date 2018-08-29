#include "mpi.h"
#include "communicator3D.h"
#include "domain3D.h"
#include "cellList3D.h"
#include "particle3D.h"
#include "rand.h"
#include "run.h"
#include "exporter.h"
using namespace std;

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  double L = 16;
  Vec_3<double> gl_l(L, L, L);
  double eta = 0.2;
  double eps = 0;
  double rho0 = 1.;
  unsigned long long seed = 12;
  Ranq1 myran(seed + my_rank);
  int n_step = 200;
  double v0 = 0.5;
  double r_cut = 1.;
  int gl_par_num = int(gl_l.x * gl_l.y * gl_l.z * rho0);

  int log_dt = 100;
  int field_dt = 100;
  int snap_dt = 100;

  typedef BiNode<VicsekPar_3> node_t;
  CellListNode_3<node_t> *cl;
  std::vector<node_t> p_arr;
  Domain_3 dm(gl_l, &cl, r_cut);
  ini_rand(p_arr, gl_par_num, myran, *cl, dm);

  // pair-wise interaction
  auto interact = [&dm, cl] (std::vector<node_t> &par_arr) {
    cal_force(par_arr, *cl, dm);
  };

  // output
  ini_output(gl_par_num, eta, eps, n_step, seed, gl_l, dm.domain_sizes());
  LogExporter *log = nullptr;
  if (my_rank == 0) {
    log = new LogExporter(log_dt);
  }
  std::cout << "succeed in initializing log for proc " << my_rank  << std::endl;

  int box_len = 1;
  if (L >= 64) {
    box_len = 2;
  }
  FieldExporter *field_exp = nullptr;
  field_exp = new FieldExporter(field_dt, 0, box_len, dm.domain_rank(),
                                dm.gl_cells_size(), dm.cells_size());
  ParticleExporter pexp(snap_dt, 0);

  OrderParaExporter order_para(100);

  auto out = [log, my_rank, field_exp, &pexp, &order_para](int i, const std::vector<node_t> &par_arr) {
    if (my_rank == 0) {
      log->record(i);
    }   
    field_exp->dump(i, par_arr);
    pexp.dump(i, par_arr);
    order_para.dump(i, par_arr);
  };
  std::cout << "succeed in initializing output for proc " << my_rank  << std::endl;
  //run(p_arr, gl_par_num, interact, integ, n_step, 1);

  if (eps == 0) {
    auto integ = [&dm, cl, &myran, eta, v0](std::vector<node_t> &par_arr) {
      integrate(par_arr, myran, *cl, dm, eta, v0);
    };
    run(p_arr, gl_par_num, interact, integ, out, n_step);
  } else {
    RandTorqueArr torques(eps, myran, dm.origin(), dm.cells_size(), dm.gl_cells_size());
    auto integ = [&dm, cl, &myran, eta, v0, &torques](std::vector<node_t> &par_arr) {
      integrate(par_arr, myran, *cl, dm, eta, v0, torques);
    };
    run(p_arr, gl_par_num, interact, integ, out, n_step);
  }

  delete log;
  delete field_exp;
  MPI_Finalize();
}
