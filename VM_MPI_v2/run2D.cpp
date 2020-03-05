#include "run2D.h"
#include "disorder2D.h"

int ini_my_par_num(int gl_par_num, MPI_Comm group_comm) {
#ifdef USE_MPI
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
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

void set_particle_num(int gl_par_num, int& my_par_num, int& my_par_num_max, MPI_Comm group_comm) {
#ifdef USE_MPI
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
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

void run_quenched(int gl_par_num, const Vec_2<double>& gl_l, double eta, double eps,
                  unsigned long long seed, int n_step, int snap_interval,
                  MPI_Comm group_comm, MPI_Comm root_comm) {
  typedef BiNode<Bird_2> node_t;
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
  Ranq1 myran(seed + my_rank);
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;
  const double v0 = 0.5;
  double rho_0 = gl_par_num / (gl_l.x * gl_l.y);

  Vec_2<int> proc_size = decompose_domain(gl_l, group_comm);
  PeriodicDomain_2 dm(gl_l, proc_size, group_comm);
  Grid_2 grid(dm, r_cut);
  CellListNode_2<node_t> cl(dm, grid);
  Communicator_2 comm(dm, grid, rho_0, 20.);

  // initialize random torques
#ifdef RANDOM_TORQUE
  RandTorque_2 disorder(eps, myran, grid, group_comm);
#elif defined RANDOM_FIELD
  RandField_2 disorder(eps, myran, grid, group_comm);
#endif

  // initialize the location and orientation of particles
  ini_rand(p_arr, gl_par_num, myran, cl, dm);

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
  auto single_move = [eta2PI, v0, &myran, &dm, &disorder, eps](node_t& p) {
    double noise = (myran.doub() - 0.5) * eta2PI;
    if (eps > 0.) {
#ifdef RANDOM_TORQUE
      noise += disorder.get_torque(p);
#elif defined RANDOM_FIELD
      disorder.apply_field(p);
#endif
    }
    move_forward(p, v0, noise, dm);
  };

  // output setting
  if (my_rank == 0) {
    mkdir("data");
#ifdef _MSC_VER
    mkdir("data\\snap");
#else
    mkdir("data/snap");
#endif
  }
  MPI_Barrier(group_comm);
  char logfile[100];
  char phifile[100];
  char snapfile[100];
 
#ifdef RANDOM_TORQUE
  snprintf(logfile, 100, "data%s%d.%d.%d.%llu.log",
    exporter::delimiter.c_str(), int(gl_l.x), int(eta * 1000), int(eps * 1000), seed);
  snprintf(phifile, 100, "data%sp%d.%d.%d.%llu.dat",
    exporter::delimiter.c_str(), int(gl_l.x), int(eta * 1000), int(eps * 1000), seed);
  snprintf(snapfile, 100, "data/snap%ss%d.%d.%d.%llu",
    exporter::delimiter.c_str(), int(gl_l.x), int(eta * 1000), int(eps * 1000), seed);
#elif defined RANDOM_FIELD
  snprintf(logfile, 100, "data%srf_%d_%g_%g_%llu.log",
    exporter::delimiter.c_str(), int(gl_l.x), eta, eps, seed);
  snprintf(phifile, 100, "data%sphi_rf_%d_%g_%g_%llu.dat",
    exporter::delimiter.c_str(), int(gl_l.x), eta, eps, seed);
  snprintf(snapfile, 100, "data/snap%ssrf_%d_%g_%g_%llu",
    exporter::delimiter.c_str(), int(gl_l.x), eta, eps, seed);
#endif
  exporter::LogExporter log(logfile, 0, n_step*2, 10000, gl_par_num, group_comm);
  exporter::OrderParaExporter_2 order_ex(phifile, 0, n_step*2, 100, gl_l, group_comm);
  exporter::SnapExporter snap_ex(snapfile, 0, n_step*2, snap_interval, group_comm);

  if (my_rank == 0) {
    log.fout << "eta=" << eta << "\n";
    log.fout << "Lx=" << gl_l.x << "\n";
    log.fout << "Ly=" << gl_l.y << "\n";
    log.fout << "particle number=" << gl_par_num << std::endl;
  }

  auto out = [&log, &order_ex, gl_par_num, &snap_ex](int i, std::vector<node_t>& par_arr) {
    log.record(i);
    order_ex.dump(i, par_arr, gl_par_num);
    snap_ex.dump(i, par_arr);
  };


  if (root_comm == MPI_COMM_WORLD) {
    for (int i = 1; i <= n_step; i++) {
      //cal_force(p_arr, cl, comm, dm);
      cal_force(p_arr, cl, comm, for_all_pair_force_fast);
      integrate(p_arr, cl, single_move, comm);
      out(i, p_arr);
    }
  } else {
    int i_step = 1;
    MPI_Win win_root;
    int finished_group = 0;
    int n_group;
    if (my_rank == 0) {
      MPI_Win_create(&i_step, 1, sizeof(int), MPI_INFO_NULL, root_comm, &win_root);
      MPI_Comm_size(root_comm, &n_group);
      std::cout << "number of groups: " << n_group << std::endl;
    }
    MPI_Bcast(&n_group, 1, MPI_INT, 0, group_comm);
    while (finished_group < n_group) {
      cal_force(p_arr, cl, comm, for_all_pair_force_fast);
      integrate(p_arr, cl, single_move, comm);
      out(i_step, p_arr);
      i_step++;
      if (i_step > n_step&& i_step % 100 == 1) {
        if (my_rank == 0) {
          finished_group = 0;
          for (int j = 0; j < n_group; j++) {
            int remote_i_step;
            MPI_Win_lock(MPI_LOCK_SHARED, j, 0, win_root);
            MPI_Get(&remote_i_step, 1, MPI_INT, j, 0, 1, MPI_INT, win_root);
            MPI_Win_unlock(j, win_root);
            if (remote_i_step > n_step) {
              finished_group++;
            }
          }
        }
        MPI_Bcast(&finished_group, 1, MPI_INT, 0, group_comm);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  std::cout << "finish" << std::endl;
}

void run_mult_bands(int gl_par_num, const Vec_2<double>& gl_l, double eta, double eps,
                    unsigned long long seed, int n_step, MPI_Comm group_comm) {
  typedef BiNode<Bird_2> node_t;
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
  Ranq1 myran(seed + my_rank);
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;
  const double v0 = 0.5;
  double rho_0 = gl_par_num / (gl_l.x * gl_l.y);

  Vec_2<int> proc_size(tot_proc, 1);
  PeriodicDomain_2 dm(gl_l, proc_size, group_comm);
  Grid_2 grid(dm, r_cut);
  CellListNode_2<node_t> cl(dm, grid);
  Communicator_2 comm(dm, grid, rho_0, 20.);

  // initialize random torques
  RandTorque_2 torque(eps, myran, grid, group_comm);

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
    exporter::delimiter.c_str(), int(eta * 1000), int(eps * 10000), gl_l.x, gl_l.y, seed);
  snprintf(phifile, 100, "data%sp%d.%d.%g.%g.%llu.dat",
    exporter::delimiter.c_str(), int(eta * 1000), int(eps * 10000), gl_l.x, gl_l.y, seed);
  snprintf(snapfile, 100, "data/snap%ss%d.%d.%g.%g.%llu",
    exporter::delimiter.c_str(), int(eta * 1000), int(eps * 10000), gl_l.x, gl_l.y, seed);
  snprintf(rhoxfile, 100, "data%srhox_%d.%d.%g.%g.%llu.bin",
    exporter::delimiter.c_str(), int(eta * 1000), int(eps * 10000), gl_l.x, gl_l.y, seed);

  exporter::LogExporter log(logfile, 0, n_step, 10000, gl_par_num, group_comm);
  exporter::OrderParaExporter_2 order_ex(phifile, 0, n_step, 100, gl_l, group_comm);
  exporter::SnapExporter snap_ex(snapfile, 0, n_step, 200000, group_comm);
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
