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

void run_birth_death(int gl_par_num, const Vec_2<double>& gl_l,
                     double eta, double alpha, unsigned long long seed,
                     int n_step, int snap_interval, int snap_block_size,
                     int box_len_birth) {
  typedef BiNode<VicsekPar_2> node_t;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  Ranq1 my_ran(seed + my_rank);
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;
  const double v0 = 0.5;
  double rho_0 = gl_par_num / (gl_l.x * gl_l.y);

  Vec_2<int> domains_size = decompose_domain(gl_l);
  Domain_2 dm(gl_l, domains_size);
  Vec_2<int> cells_size{};
  Vec_2<double> cell_len{};
  CellListBase_2::partition(dm.l(), r_cut, cells_size, cell_len);
  Vec_2<int> gl_cells_size(cells_size.x * domains_size.x, cells_size.y * domains_size.y);
  CellListNode_2<node_t> cl(cells_size, cell_len, gl_l, dm.origin(), dm.flag_comm());
  Communicator comm(gl_l, rho_0, 20., cells_size, domains_size);
  
  ini_rand(p_arr, gl_par_num, my_ran, cl, dm);

  std::vector<unsigned short> n_arr(cl.n_cells(), 0);
  std::vector<Vec_2<double>> v_arr(cl.n_cells(), Vec_2<double>());

  // pair-wise interaction
  auto interact = [&comm, &cl](std::vector<node_t> &par_arr) {
    cal_force(par_arr, cl, comm);
  };

#ifdef OUTPUT_ON
  ini_output(gl_par_num, eta, alpha, n_step, seed, gl_l, domains_size);
  LogExporter *log_ex = nullptr;
  if (my_rank == 0) {
    log_ex = new LogExporter(10000);
  }
  OrderParaExporter phi_ex(100);
  const int field_record_dt = snap_interval;
  const int cg_box_size = snap_block_size;
  int n_frames_per_file = 20;
  int count_frames = 0;
  int n_steps_per_file = n_frames_per_file * field_record_dt;
  FieldExporter *field_ex = new FieldExporter(field_record_dt, cg_box_size, dm.rank(),
                                              gl_cells_size, cells_size);
  auto out = [log_ex, &phi_ex, &field_ex, my_rank, n_frames_per_file,
              &count_frames, field_record_dt, cg_box_size, &dm, &cl, &n_arr,
              &gl_cells_size, &cells_size](int i, std::vector<node_t> &par_arr) {
    if (my_rank == 0) {
      log_ex->record(i);
    }
    phi_ex.dump(i, par_arr, cl, n_arr);
    if (field_ex->dump(i, par_arr)) {
      count_frames++;
      if (count_frames % n_frames_per_file == 0) {
        delete field_ex;
        field_ex = new FieldExporter(field_record_dt,
          cg_box_size, dm.rank(), gl_cells_size, cells_size);
      }
    }
  };
#endif
  auto single_move = [eta, v0, &my_ran, &dm](node_t &p) {
    p.move(eta, v0, my_ran);
  };
  auto f_rate = [rho_0, alpha](double rho) {
    return (1. - rho / rho_0) * alpha;
  };

  Vec_2<int> box_birth(box_len_birth, box_len_birth);
  auto integ = [&comm, &cl, single_move, &n_arr, &v_arr, f_rate,
                &my_ran, &box_birth](std::vector<node_t> &par_arr) {
    integrate(par_arr, cl, single_move, comm, n_arr, v_arr);
    cl.birth_death(par_arr, f_rate, my_ran, n_arr, v_arr, box_birth);
  };

  run(p_arr, interact, integ, out, n_step);


#ifdef OUTPUT_ON
  delete log_ex;
#endif
}

void run_test(int gl_par_num, const Vec_2<double>& gl_l, double eta, unsigned long long seed, int n_step) {
  typedef BiNode<VicsekPar_2> node_t;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  Ranq1 my_ran(seed + my_rank);
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;
  const double v0 = 0.5;
  double rho_0 = gl_par_num / (gl_l.x * gl_l.y);

  Vec_2<int> domains_size = decompose_domain(gl_l);
  Domain_2 dm(gl_l, domains_size);
  Vec_2<int> cells_size{};
  Vec_2<double> cell_len{};
  CellListBase_2::partition(dm.l(), r_cut, cells_size, cell_len);
  Vec_2<int> gl_cells_size(cells_size.x * domains_size.x, cells_size.y * domains_size.y);
  CellListNode_2<node_t> cl(cells_size, cell_len, gl_l, dm.origin(), dm.flag_comm());
  Communicator comm(gl_l, rho_0, 20., cells_size, domains_size);

  ini_rand(p_arr, gl_par_num, my_ran, cl, dm);

  auto single_move = [eta, v0, &my_ran, &dm](node_t &p) {
    p.move(eta, v0, my_ran, dm);
  };

  for (int i = 1; i <= n_step; i++) {
    cal_force(p_arr, cl, comm, dm);
    integrate(p_arr, cl, single_move, comm);
    if (i % 1000 == 0) {
      double vel_mean[2];
      get_mean_vel(vel_mean, p_arr, gl_par_num, true);
      if (my_rank == 0) {
        double phi = std::sqrt(vel_mean[0] * vel_mean[0] + vel_mean[1] * vel_mean[1]);
        std::cout << "phi = " << phi << std::endl;
      }
    }
  }
  std::cout << "finish" << std::endl;
}
