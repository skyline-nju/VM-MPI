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

void run_rand_torque(int gl_par_num, const Vec_2<double> &gl_l, double eta,
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
    log_ex = new LogExporter(1000);
  }
  OrderParaExporter phi_ex(100);
  const int field_record_dt = 1000;
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

  if (eps == 0) {
    auto single_move = [eta, v0, &myran](node_t &p, const Domain_2& domain) {
      p.move(eta, v0, myran, domain);
    };
    auto integ = [&dm, cl, single_move](std::vector<node_t> &par_arr) {
      integrate(par_arr, *cl, dm, single_move);
    };
#ifdef OUTPUT_ON
    run(p_arr, gl_par_num, interact, integ, out, n_step);
#else
    run(p_arr, gl_par_num, interact, integ, n_step, 100);
#endif
  } else {
    RandTorqueArr_2 torques(eps, myran, dm.origin(), dm.cells_size(), dm.gl_cells_size());
    auto single_move = [eta, v0, &myran, &torques](node_t &p, const Domain_2& domain) {
      double tau = torques.get_torque(p);
      p.move_torque(eta, v0, tau, myran, domain);
    };
    auto integ = [&dm, cl, single_move](std::vector<node_t> &par_arr) {
      integrate(par_arr, *cl, dm, single_move);
    };
#ifdef OUTPUT_ON
    run(p_arr, gl_par_num, interact, integ, out, n_step);
#else
    run(p_arr, gl_par_num, interact, integ, n_step, 100);
#endif
  }

#ifdef OUTPUT_ON
  delete log_ex;
#endif
}

void run_osc_field(int gl_par_num, const Vec_2<double>& gl_l, double eta,
                   double h0, int t_half, int n_period, unsigned long long seed,
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

  auto single_move = [eta, v0, &myran](node_t &p, const Domain_2& domain, double h) {
    p.move_field(eta, v0, h, myran, domain);
  };

  int t = 0;
  int t_equil = 20000;
  int t_period = 2 * t_half;
  auto integ = [&dm, cl, single_move, &t, t_equil, t_half, t_period, h0](std::vector<node_t> &par_arr) {
    double h = h0;
    if (t >= t_equil && (t - t_equil) % t_period < t_half) {
      h = -h0;
    }
    integrate(par_arr, *cl, dm, single_move, h);
    t++;
  };

  const int n_steps = t_equil + t_half * 2 * n_period;
#ifdef OUTPUT_ON
  ini_output(gl_par_num, eta, h0, t_half, t_equil, n_steps, seed, gl_l, dm.domain_sizes());
  LogExporter *log_ex = nullptr;
  if (my_rank == 0) {
    log_ex = new LogExporter(1000);
  }
  OrderParaExporter phi_ex(1);
  const int field_record_dt = 1000;
  const int cg_box_size = 1;
  FieldExporter field_ex(field_record_dt, 0, cg_box_size,
    dm.domain_rank(), dm.gl_cells_size(), dm.cells_size());
  auto out = [log_ex, &phi_ex, &field_ex, my_rank](int i, std::vector<node_t> &par_arr) {
    if (my_rank == 0) {
      log_ex->record(i);
    }
    phi_ex.dump(i, par_arr);
    field_ex.dump(i, par_arr);
  };
  run(p_arr, gl_par_num, interact, integ, out, n_steps);
#else
  run(p_arr, gl_par_num, interact, integ, n_steps, 100);
#endif
}
