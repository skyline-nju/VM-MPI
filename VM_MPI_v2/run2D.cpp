#include "run2D.h"
#include "disorder2D.h"

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

void run_mult_bands(int gl_par_num, const Vec_2<double>& gl_l, double eta, double eps, unsigned long long seed, int n_step) {
  typedef BiNode<Bird_2> node_t;
  int my_rank, tot_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  Ranq1 myran(seed + my_rank);
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;
  const double v0 = 0.5;
  double rho_0 = gl_par_num / (gl_l.x * gl_l.y);

  Vec_2<int> proc_size(tot_proc, 1);
  PeriodicDomain_2 dm(gl_l, proc_size);
  Grid_2 grid(dm, r_cut);
  CellListNode_2<node_t> cl(dm, grid);
  Communicator_2 comm(dm, grid, rho_0, 20.);

  // initialize random torques, the size of each subdomain must be equal
  RandTorque_2 torque(eps, myran, dm.origin(), grid.n(), grid.gl_n());

  // initialize the location and orientation of particles
  ini_rand(p_arr, gl_par_num, myran, cl, dm);
  Vec_2<double> u0(1., 0);
  for (auto& i : p_arr) {
    i.ori = i.ori_next = u0;
  }

  // cal force
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

  // integrate
  double eta2PI = eta * 2.0 * PI;
  auto single_move = [eta2PI, v0, &myran, &dm, &torque, eps](node_t& p) {
    double noise = (myran.doub() - 0.5) * eta2PI;
    // random torque
    if (eps > 0.) {
      noise += torque.get_torque(p);
    }
    move_forward(p, v0, noise, dm);
  };

  // output setting
  mkdir("data");
  mkdir("data/snap");
  char logfile[100];
  char phifile[100];
  char snapfile[100];
  char rhoxfile[100];

  snprintf(logfile, 100, "data%s%d.%d.%g.%g.%llu.log",
    exporter::delimiter.c_str(), int(eta * 1000), int(eps * 1000), gl_l.x, gl_l.y, seed);
  snprintf(phifile, 100, "data%sp%d.%d.%g.%g.%llu.dat",
    exporter::delimiter.c_str(), int(eta * 1000), int(eps * 1000), gl_l.x, gl_l.y, seed);
  snprintf(snapfile, 100, "data/snap%ss%d.%d.%g.%g.%llu",
    exporter::delimiter.c_str(), int(eta * 1000), int(eps * 1000), gl_l.x, gl_l.y, seed);
  snprintf(rhoxfile, 100, "data%srhox_%d.%d.%g.%g.%llu.bin",
    exporter::delimiter.c_str(), int(eta * 1000), int(eps * 1000), gl_l.x, gl_l.y, seed);

  exporter::LogExporter log(logfile, 0, n_step, 10000, gl_par_num);
  exporter::OrderParaExporter_2 order_ex(phifile, 0, n_step, 100);
  exporter::SnapExporter snap_ex(snapfile, 0, n_step, 200000);
  exporter::RhoxExporter rhox_ex(rhoxfile, 0, n_step, 100, grid, dm);

  if (my_rank == 0) {
    log.fout << "eta=" << eta << "\n";
    log.fout << "Lx=" << gl_l.x << "\n";
    log.fout << "Ly=" << gl_l.y << "\n";
    log.fout << "particle number=" << gl_par_num << std::endl;
  }

  auto out = [&log, &order_ex, gl_par_num, &snap_ex, &rhox_ex](int i, std::vector<node_t>& par_arr) {
    log.record(i);
    order_ex.dump(i, par_arr, gl_par_num);
    snap_ex.dump(i, par_arr);
    rhox_ex.dump(i, par_arr);
  };

  // run
  for (int i = 1; i <= n_step; i++) {
    //cal_force(p_arr, cl, comm, dm);
    cal_force(p_arr, cl, comm, for_all_pair_force_fast);
    integrate(p_arr, cl, single_move, comm);
    out(i, p_arr);
  }
  std::cout << "finish" << std::endl;

}
