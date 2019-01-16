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

//void run(int gl_par_num, const Vec_2<double>& gl_l,
//         double eta, double eps, unsigned long long seed,
//         int n_step, int n_save_snap, int block_size) {
//  typedef BiNode<VicsekPar_2> node_t;
//  int my_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//  Ranq1 myran(seed + my_rank);
//  CellListNode_2<node_t> *cl;
//  std::vector<node_t> p_arr;
//  const double r_cut = 1.;
//  const double v0 = 0.5;
//  Domain_2 dm(gl_l, &cl, r_cut);
//  ini_rand(p_arr, gl_par_num, myran, *cl, dm);
//
//  // pair-wise interaction
//  auto interact = [&dm, cl](std::vector<node_t> &par_arr) {
//    cal_force(par_arr, *cl, dm);
//  };
//
//#ifdef OUTPUT_ON
//  ini_output(gl_par_num, eta, eps, n_step, seed, gl_l, dm.domain_sizes());
//  LogExporter *log_ex = nullptr;
//  if (my_rank == 0) {
//    log_ex = new LogExporter(10000);
//  }
//  OrderParaExporter phi_ex(100);
//  const int field_record_dt = n_save_snap;
//  const int cg_box_size = block_size;
//  FieldExporter field_ex(field_record_dt,  cg_box_size,
//                         dm.domain_rank(), dm.gl_cells_size(), dm.cells_size());
//  auto out = [log_ex, &phi_ex, &field_ex, my_rank](int i, std::vector<node_t> &par_arr) {
//    if (my_rank == 0) {
//      log_ex->record(i);
//    }
//    phi_ex.dump(i, par_arr);
//    field_ex.dump(i, par_arr);
//  };
//#endif
//#ifdef DISORDER_ON
//  Disorder_2 disorder(eps, myran, dm.origin(), dm.cells_size(), dm.gl_cells_size());
//  auto single_move = [eta, v0, &myran, &disorder](node_t &p, const Domain_2& domain) {
//    p.move(eta, v0, myran, domain, disorder);
//  };
//#else
//  auto single_move = [eta, v0, &myran](node_t &p, const Domain_2& domain) {
//    p.move(eta, v0, myran, domain);
//  };
//#endif
//  auto integ = [&dm, cl, single_move](std::vector<node_t> &par_arr) {
//    integrate(par_arr, *cl, dm, single_move);
//  };
//
//#ifdef OUTPUT_ON
//  run(p_arr, gl_par_num, interact, integ, out, n_step);
//#else
//  run(p_arr, gl_par_num, interact, integ, n_step, 100);
//#endif
//
//#ifdef OUTPUT_ON
//  delete log_ex;
//#endif
//}

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

  std::vector<unsigned short> n_arr(cl.ncells(), 0);
  std::vector<Vec_2<double>> v_arr(cl.ncells(), Vec_2<double>());

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
    p.move(eta, v0, my_ran, dm);
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
