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
    my_par_num_max = 15 * my_par_num;
  } else {
    my_par_num_max = 1.1 * my_par_num;
  }
#else
  my_par_num = gl_par_num;
  my_par_num_max = gl_par_num;
#endif
}

void run(int gl_par_num, const Vec_2<double>& gl_l,
         double eta, double eps, double v0,
         int n_step, std::string& ini_mode,
         unsigned long long seed, int seed2,
         int snap_dt, int field_dt, int field_dx,
         MPI_Comm group_comm, MPI_Comm root_comm) {
  typedef BiNode<Bird_2> node_t;
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
  Ranq1 myran(seed * 10 + my_rank);
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;
  double rho_0 = gl_par_num / (gl_l.x * gl_l.y);
  double eta2PI = eta * 2.0 * PI;

  //Vec_2<int> proc_size = decompose_domain(gl_l, group_comm);
  Vec_2<int> proc_size = Vec_2<int>(1, 1);
  PeriodicDomain_2 dm(gl_l, proc_size, group_comm);
  Grid_2 grid(dm, r_cut);
  CellListNode_2<node_t> cl(dm, grid);
  Communicator_2 comm(dm, grid, rho_0, 80.);

  // initialize aligning force
  auto f1 = [](node_t* p1, node_t* p2) {
    polar_align(*p1, *p2, p2->pos - p1->pos);
};
  auto f2 = [&dm](node_t* p1, node_t* p2) {
    polar_align(*p1, *p2, dm);
  };
  auto f3 = [](node_t* p1, node_t* p2, const Vec_2<double>& offset) {
    polar_align(*p1, *p2, p2->pos - p1->pos + offset);
  };
#if defined REF_WALL_Y || defined REF_WALL_XY
  auto for_all_pair_force = [&cl, &f1, &f2]() {
    cl.for_each_pair_slow(f1, f2);
  };
#else
  auto for_all_pair_force = [&cl, &f1, &f3]() {
    cl.for_each_pair_fast(f1, f3);
  };
#endif

  // initialize particles
  int t_beg;
  ini_particles(gl_par_num, p_arr, ini_mode, myran, seed2, cl, dm, t_beg);

  // initialize integrator
  Ranq1 myran2(seed + t_beg + my_rank);
  auto single_move = [eta2PI, v0, &myran2, &dm, eps](node_t& p) {
    double noise = (myran2.doub() - 0.5) * eta2PI;
    move_forward(p, v0, noise, dm);
  };
  
  // output setting
  char folder[255];
  char snapfolder[255];
#ifndef _MSC_VER
  snprintf(folder, 255, "/scratch03.local/yduan/polar_LG_W/%g_%g/", gl_l.x, gl_l.y);
  snprintf(snapfolder, 255, "%ssnap/", folder);
#else
  snprintf(folder, 255, "data\\%g_%g\\", gl_l.x, gl_l.y);
  snprintf(snapfolder, 255, "%ssnap\\", folder);
#endif
  if (my_rank == 0) {
    mkdir(folder);
    mkdir(snapfolder);
  }
  MPI_Barrier(group_comm);
  char logfile[255];
  char phifile[255];
  char snapfile[255];
  char fieldfile[255];
  char basename[255];

  if (gl_l.x == gl_l.y) {
    snprintf(basename, 255, "%g_%.3f_%.3f_%.1f_%llu", gl_l.x, eta, rho_0, v0, seed);
  } else {
    snprintf(basename, 255, "%g_%g_%.3f_%.3f_%.1f_%llu", gl_l.x, gl_l.y, eta, rho_0, v0, seed);
  }

  snprintf(logfile, 255, "%s%s_%d.log", folder, basename, t_beg);
  snprintf(phifile, 255, "%s%s_%d.dat", folder, basename, t_beg);
  snprintf(snapfile, 255, "%ss%s", snapfolder, basename);
  snprintf(fieldfile, 255, "%s%s_%d_%d_%d.bin", folder, basename, field_dx, field_dt, t_beg);

  exporter::LogExporter log(logfile, 0, n_step * 2, 10000, gl_par_num, group_comm);
  exporter::OrderParaExporter_2 order_ex(phifile, 0, n_step * 2, 100, gl_l, group_comm);
  exporter::SnapExporter snap_ex(snapfile, 0, n_step * 2, snap_dt, t_beg, group_comm);
  exporter::FeildExporter field_ex(fieldfile, 0, n_step * 2, field_dt, grid, dm, field_dx);

  if (my_rank == 0) {
    log.fout << "eta=" << eta << "\n";
    log.fout << "rho_0=" << rho_0 << "\n";
    log.fout << "v_0=" << v0 << "\n";
    log.fout << "Lx=" << gl_l.x << "\n";
    log.fout << "Ly=" << gl_l.y << "\n";
    log.fout << "particle number=" << gl_par_num << "\n";
    log.fout << "core number=" << tot_proc << "\n";
    log.fout << "field dx=" << field_dx << "\n";
    log.fout << "field dt=" << field_dt << "\n";
    log.fout << "t_beg=" << t_beg << std::endl;
  }

  auto out = [&log, &order_ex, gl_par_num, &snap_ex, &field_ex](int i, std::vector<node_t>& par_arr) {
    log.record(i);
    order_ex.dump(i, par_arr, gl_par_num);
    snap_ex.dump(i, par_arr);
    field_ex.dump(i, par_arr);
  };

  // run

  bool thick_shell = v0 > 1.;
  int gl_tot_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &gl_tot_proc);
  if (gl_tot_proc == tot_proc) {
    for (int i = 1; i <= n_step; i++) {
      cal_force(p_arr, cl, comm, for_all_pair_force);
      integrate(p_arr, cl, single_move, comm, thick_shell);

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
      cal_force(p_arr, cl, comm, for_all_pair_force);
      integrate(p_arr, cl, single_move, comm, thick_shell);
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
