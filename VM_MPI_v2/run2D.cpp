#include "run2D.h"

int ini_my_par_num(int gl_par_num) {
#ifdef USE_MPI
  int my_rank, tot_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  int my_par_num;
  if (my_rank > 0) {
    my_par_num = gl_par_num / tot_proc;
  } else {
    my_par_num = gl_par_num - gl_par_num / tot_proc * (tot_proc - 1);
  }
  return my_par_num;
#else
  return gl_par_num;
#endif
}

void run(int gl_par_num, const Vec_2<double>& gl_l,
         double eta, double eps, unsigned long long seed,
         int n_step, int n_save_snap, int block_size,
         double r_cut, double v0) {
  typedef BiNode<VicsekPar_2> node_t;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  Ranq1 myran(seed + my_rank);
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
  LogExporter *log_ex = nullptr;
  if (my_rank == 0) {
    log_ex = new LogExporter(10000);
  }
  OrderParaExporter phi_ex(100);
  const int field_record_dt = n_save_snap;
  const int cg_box_size = block_size;
  FieldExporter field_ex(field_record_dt, 0, cg_box_size,
                         dm.domain_rank(), dm.gl_cells_size(), dm.cells_size());
  auto out = [log_ex, &phi_ex, &field_ex, my_rank](int i, std::vector<node_t> &par_arr) {
    if (my_rank == 0) {
      log_ex->record(i);
    }
    phi_ex.dump(i, par_arr);
    field_ex.dump(i, par_arr);
  };
#endif
#ifdef DISORDER_ON
  Disorder_2 disorder(eps, myran, dm.origin(), dm.cells_size(), dm.gl_cells_size());
  auto single_move = [eta, v0, &myran, &disorder](node_t &p, const Domain_2& domain) {
    p.move(eta, v0, myran, domain, disorder);
  };
#else
  auto single_move = [eta, v0, &myran, &disorder](node_t &p, const Domain_2& domain) {
    p.move(eta, v0, myran, domain);
  };
#endif
  auto integ = [&dm, cl, single_move](std::vector<node_t> &par_arr) {
    integrate(par_arr, *cl, dm, single_move);
  };

#ifdef OUTPUT_ON
  run(p_arr, gl_par_num, interact, integ, out, n_step);
#else
  run(p_arr, gl_par_num, interact, integ, n_step, 100);
#endif

#ifdef OUTPUT_ON
  delete log_ex;
#endif
}

#ifdef DENSITY_NOISE
void run_density_noise(int gl_par_num, const Vec_2<double> &gl_l, double eta,
                       double eps, int n_step, unsigned long long seed,
                       double r_cut, double v0) {
  typedef BiNode<VicsekPar_2> node_t;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  Ranq1 myran(seed + my_rank);
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
  LogExporter *log_ex = nullptr;
  if (my_rank == 0) {
    log_ex = new LogExporter(10000);
  }
  OrderParaExporter phi_ex(100);
  const int field_record_dt = 50000;
  const int cg_box_size = 2;
  FieldExporter field_ex(field_record_dt, 0, cg_box_size,
                         dm.domain_rank(), dm.gl_cells_size(), dm.cells_size());
  auto out = [log_ex, &phi_ex, &field_ex, my_rank](int i, std::vector<node_t> &par_arr) {
    if (my_rank == 0) {
      log_ex->record(i);
    }
    phi_ex.dump(i, par_arr);
    field_ex.dump(i, par_arr);
  };
#endif

  auto single_move = [eta, v0, &myran, eps](node_t &p, const Domain_2& domain) {
    p.move_density_noise(eta, v0, eps, myran, domain);
  };
  auto integ = [&dm, cl, single_move](std::vector<node_t> &par_arr) {
    integrate(par_arr, *cl, dm, single_move);
  };
#ifdef OUTPUT_ON
  run(p_arr, gl_par_num, interact, integ, out, n_step);
#else
  run(p_arr, gl_par_num, interact, integ, n_step, 100);
#endif
  

#ifdef OUTPUT_ON
  delete log_ex;
#endif
}
#endif