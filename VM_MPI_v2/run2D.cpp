#include "run2D.h"
#include "exporter2D_b.h"

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

void set_particle_num(int gl_par_num, int& my_par_num, int& my_par_num_max) {
#ifdef USE_MPI
  int my_rank, tot_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  if (my_rank > 0) {
    my_par_num = gl_par_num / tot_proc;
  } else {
    my_par_num = gl_par_num - gl_par_num / tot_proc * (tot_proc - 1);
  }
  my_par_num = my_par_num;
  if (tot_proc > 1) {
    my_par_num_max = 3 * my_par_num;
  } else {
    my_par_num_max = 1.1 * my_par_num;
  }
#else
  my_par_num = gl_par_num;
  my_par_num_max = gl_par_num;
#endif
}

void run(int gl_par_num, const Vec_2<double>& gl_l, double eta,
         unsigned long long seed, int n_step, const char* file_in) {
  typedef BiNode<Bird_2> node_t;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  Ranq1 myran(seed + my_rank);
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;
  const double v0 = 0.5;
  double rho_0 = gl_par_num / (gl_l.x * gl_l.y);

  int tot_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  Vec_2<int> proc_size(tot_proc, 1);
  //SubDomain_2 dm(gl_l, proc_size, r_cut);
  PeriodicDomain_2 dm(gl_l, proc_size);
  Grid_2 grid(dm, r_cut);

  CellListNode_2<node_t> cl(dm, grid);
  //Communicator_2 comm(dm, rho_0, 20.);
  Communicator_2 comm(dm, grid, rho_0, 20.);

  if (file_in) {
    ini_from_snap(p_arr, gl_par_num, cl, dm, file_in);
  } else {
    ini_rand(p_arr, gl_par_num, myran, cl, dm);
  }
  
  std::cout << p_arr.size() << " particles at proc " << my_rank << std::endl;
       
  auto f1 = [](node_t* p1, node_t* p2) {
    polar_align(*p1, *p2, p2->pos - p1->pos);
  };

  auto f2 = [&dm](node_t* p1, node_t* p2) {
    polar_align(*p1, *p2, dm);
  };

  auto f3 = [](node_t* p1, node_t* p2, const Vec_2<double>& offset) {
    polar_align(*p1, *p2, p2->pos - p1->pos + offset);
  };

  auto for_all_pair_force_slow = [&cl, &f1, &f2]() {
    cl.for_each_pair_slow(f1, f2);
  };

  auto for_all_pair_force_fast = [&cl, &f1, &f3]() {
    cl.for_each_pair_fast(f1, f3);
  };

  double eta2PI = eta * 2.0 * PI;
  auto single_move = [eta2PI, v0, &myran, &dm](node_t &p) {
    double noise = (myran.doub() - 0.5) * eta2PI;
    move_forward(p, v0, noise, dm);
  };

  exporter_ini(gl_par_num, eta, 0., n_step, seed, gl_l, dm.proc_size());

  LogExporter *log_ex = nullptr;
  OrderParaExporter phi_ex(100);
  SnapExporter snap_ex(10000, n_step);
  RhoxExpoerter rhox_ex(10000, n_step, grid, dm.origin());

  if (my_rank == 0) {
    log_ex = new LogExporter(10000);
  }
  auto out = [log_ex, &phi_ex, my_rank, &dm, &snap_ex, &rhox_ex](int i, std::vector<node_t> &par_arr) {
    if (my_rank == 0) {
      log_ex->record(i);
    }
    phi_ex.dump(i, par_arr);
    snap_ex.dump(i, par_arr);
    rhox_ex.dump(i, par_arr);
  };

  for (int i = 1; i <= n_step; i++) {
    //cal_force(p_arr, cl, comm, dm);
    cal_force(p_arr, cl, comm, for_all_pair_force_fast);
    integrate(p_arr, cl, single_move, comm);

    out(i, p_arr);
  }
  std::cout << "finish" << std::endl;
}
