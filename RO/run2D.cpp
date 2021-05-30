#include "run2D.h"
#include "disorder2D.h"

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
    my_par_num_max = 10 * my_par_num;
  } else {
    my_par_num_max = 1.1 * my_par_num;
  }
#else
  my_par_num = gl_par_num;
  my_par_num_max = gl_par_num;
#endif
}

void run_RO(int gl_par_num, const Vec_2<double>& gl_l,
            double eta, double rho_s, double eps,
            unsigned long long seed, unsigned long long seed2,
            int n_step, std::string& ini_mode,
            MPI_Comm group_comm, MPI_Comm root_comm) {
  typedef BiNode<Bird_2> node_t;
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;
  const double v0 = 0.5;
  double rho_0 = gl_par_num / (gl_l.x * gl_l.y);

  Vec_2<int> proc_size = decompose_domain(gl_l, group_comm);
  PeriodicDomain_2 dm(gl_l, proc_size, group_comm);
  Grid_2 grid(dm, r_cut);
  CellListNode_2<node_t> cl(dm, grid);
  Communicator_2 comm(dm, grid, rho_0, 40.);

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
  char time_ave_f1[100];
  char field_f1[100];
  char scatter_f[100];
  char basename[100];
  if (gl_l.x == gl_l.y) {
    snprintf(basename, 100, "%g_%.3f_%.5f_%g_%llu_%03llu", gl_l.x, eta, rho_s, eps, seed, seed2);
  } else {
    snprintf(basename, 100, "%g_%g_%.3f_%.5f_%g_%llu_%03llu", gl_l.x, gl_l.y, eta, rho_s, eps, seed, seed2);
  }
  snprintf(scatter_f, 100, "data/scatter_%s.bin", basename);

  int ave_t_beg1 = 400000;
  int ave_dt1 = 100000;
  int field_t_beg1 = 0;
  int field_dt1 = 10000;

  // initialize random obstacles
  int n_ob = round(rho_s * gl_l.x * gl_l.y);
  Ranq1 myran(seed + my_rank);
  RandScatter disorder(dm, grid, n_ob, myran, group_comm, scatter_f);

  // initialize the location and orientation of particles
  int t_beg;
  ini_particles(gl_par_num, p_arr, ini_mode, myran, seed2, cl, dm, t_beg);

  snprintf(logfile, 100, "data/RO_%s_%d.log", basename, t_beg);
  snprintf(phifile, 100, "data/RO_%s_%d.dat", basename, t_beg);
  snprintf(snapfile, 100, "data/snap/RO_%s", basename);
  snprintf(field_f1, 100, "data/RO_field_%s_%d_%d.bin", basename, field_dt1, t_beg);
  snprintf(time_ave_f1, 100, "data/RO_ave_%s_%d_%d.bin", basename, ave_dt1, t_beg);

  exporter::LogExporter log(logfile, 0, n_step * 2, 10000, gl_par_num, group_comm);
  exporter::OrderParaExporter_2 order_ex(phifile, 0, n_step * 2, 100, gl_l, group_comm);
  exporter::SnapExporter snap_ex(snapfile, 0, n_step * 2, 10000, t_beg,  group_comm);
  exporter::FeildExporter field_ex1(field_f1, field_t_beg1, n_step * 2, field_dt1, grid, dm);
  exporter::TimeAveFeildExporter ave_ex1(time_ave_f1, ave_t_beg1, n_step * 2, ave_dt1, grid, dm);

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
  Ranq2 myran2(myran.int64() + seed2 + t_beg);
  auto single_move = [eta2PI, v0, &myran2, &dm, &disorder, eps](node_t& p) {
    disorder.scattering(p, eps);
    double noise = (myran2.doub() - 0.5) * eta2PI;
    move_forward(p, v0, noise, dm);
  };

  if (my_rank == 0) {
    log.fout << "eta=" << eta << "\n";
    log.fout << "Lx=" << gl_l.x << "\n";
    log.fout << "Ly=" << gl_l.y << "\n";
    log.fout << "particle number=" << gl_par_num << std::endl;
  }

  auto out = [&log, &order_ex, gl_par_num, &snap_ex, &ave_ex1, &field_ex1](int i, std::vector<node_t>& par_arr) {
    log.record(i);
    order_ex.dump(i, par_arr, gl_par_num);
    snap_ex.dump(i, par_arr);
    ave_ex1.dump(i, par_arr);
    field_ex1.dump(i, par_arr);
  };

  if (root_comm == MPI_COMM_WORLD) {
    for (int i = 1; i <= n_step; i++) {
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
      if (i_step > n_step && i_step % 100 == 1) {
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
