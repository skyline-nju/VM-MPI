#include "run2D.h"

#ifdef CONTINUE_DYNAMIC
#include "chiral_vicsek.h"

void run(int gl_par_num, const Vec_2<double>& gl_l,
         double D_0, double omega0,
         unsigned long long seed,
         int n_step, std::string& ini_mode,
         MPI_Comm group_comm, MPI_Comm root_comm) {
  typedef BiNode<Bird_2> node_t;
  int my_rank = 0;
  int tot_proc = 1;
#ifdef USE_MPI
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
#endif
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;
  const double v0 = 1.0;
  const double h = 0.1;
  double rho_0 = gl_par_num / (gl_l.x * gl_l.y);
  double sqrt_24Dh = std::sqrt(24 * D_0 * h);
  double v0h = v0 * h;

  Vec_2<int> proc_size = decompose_domain(gl_l, group_comm);
  PeriodicDomain_2 dm(gl_l, proc_size, group_comm);
  Grid_2 grid(dm, r_cut);
  CellListNode_2<node_t> cl(dm, grid);
#ifdef USE_MPI
  Communicator_2 comm(dm, grid, rho_0, 40.);
#endif

  exporter::create_folders(my_rank, group_comm);

  char logfile[200];
  char phifile[200];
  char snapfile[200];
  char time_ave_f1[200];
  char field_f1[200];
  char scatter_f[200];
  char basename[100];
  if (gl_l.x == gl_l.y) {
    snprintf(basename, 100, "%g_%.4f_%.4f_%.2f_%llu", gl_l.x, D_0, omega0, rho_0, seed);
  } else {
    snprintf(basename, 100, "%g_%g_%.4f_%.4f_%.2f_%llu", gl_l.x, gl_l.y, D_0, omega0, rho_0, seed);
  }
  snprintf(scatter_f, 200, "data/scatter_%s.bin", basename);


  // initialize random obstacles
  Ranq1 myran(seed + my_rank);

  // initialize the location and orientation of particles
  int t_beg;
  ini_particles(gl_par_num, p_arr, ini_mode, myran, cl, dm, t_beg);

  int field_t_beg1 = 0;
  int field_dt1 = 10000;
  int field_dx = 2;
  snprintf(logfile, 200, "data/RO_%s_%d.log", basename, t_beg);
  snprintf(phifile, 200, "data/RO_%s_%d.dat", basename, t_beg);
  snprintf(snapfile, 200, "data/snap/RO_%s", basename);
  snprintf(field_f1, 200, "data/RO_field_%s_%d_%d_%d.bin", basename, field_dx, field_dt1, t_beg);

  exporter::LogExporter log(logfile, 0, n_step * 2, 10000, gl_par_num, group_comm);
  exporter::OrderParaExporter_2 order_ex(phifile, 0, n_step * 2, 100, gl_l, group_comm);
  exporter::SnapExporter snap_ex(snapfile, 0, n_step * 2, 10000, t_beg, group_comm);
  exporter::FeildExporter field_ex1(field_f1, field_t_beg1, n_step * 2, field_dt1, grid, dm, field_dx);
 
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
  Ranq2 myran2(myran.int64() + t_beg);

  double omega0_h = omega0 * h;
  auto single_move = [h, omega0_h, sqrt_24Dh, v0h, &myran2, &dm](node_t& p) {
    double d_theta = (myran2.doub() - 0.5) * sqrt_24Dh + omega0_h;
    if (p.n_neighb > 1) {
      d_theta += p.tau / (p.n_neighb - 1) * h;
    }
    move_forward(p, v0h, d_theta, dm);
    p.tau = 0.;
    p.n_neighb = 1;
  };

  if (my_rank == 0) {
    log.fout << "D_0=" << D_0 << "\n";
    log.fout << "omega0=" << omega0 << "\n";
    log.fout << "rho_0=" << rho_0 << "\n";
    log.fout << "v0=" << v0 << "\n";
    log.fout << "dt=" << h << "\n";
    log.fout << "Lx=" << gl_l.x << "\n";
    log.fout << "Ly=" << gl_l.y << "\n";
    log.fout << "particle number=" << gl_par_num << "\n";
    log.fout << "tot_proc=" << tot_proc << std::endl;
  }

  auto out = [&log, &order_ex, gl_par_num, &snap_ex, &field_ex1](int i, std::vector<node_t>& par_arr) {
    log.record(i);
    order_ex.dump(i, par_arr, gl_par_num);
    snap_ex.dump(i, par_arr);
    field_ex1.dump(i, par_arr);
  };

#ifndef USE_MPI
  auto run_one_step = [&p_arr, &cl, &for_all_pair_force_fast, &single_move, &out](int i) {
    for_all_pair_force_fast();
    integrate2(p_arr, cl, single_move);
    out(i, p_arr);
  };
#else
  auto run_one_step = [&p_arr, &cl, &comm, &for_all_pair_force_fast, &single_move, &out](int i) {
    cal_force(p_arr, cl, comm, for_all_pair_force_fast);
    integrate2(p_arr, cl, single_move, comm);
    out(i, p_arr);
    if (i % 1000 == 0) {
      for (auto& p : p_arr) {
        p.ori.normalize();
      }
    }
  };
#endif
  run(run_one_step, n_step, group_comm, root_comm);
}


int main(int argc, char* argv[]) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  int my_rank, tot_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);

  int arg_size = 7;
  if ((argc - 1) % arg_size != 0) {
    std::cout << "Error, argc = " << argc << std::endl;
    exit(1);
  }
  int n_group = 1;
  int cores_per_group = tot_proc / n_group;
  MPI_Group group, gl_group, root_group;
  MPI_Comm group_comm, root_comm;
  int* ranks = new int[cores_per_group];  // rank of processer belong to the same group
  int* root_ranks = new int[n_group];     // rank of root processer of each group
  int my_group = my_rank / cores_per_group;
  for (int i = 0; i < cores_per_group; i++) {
    ranks[i] = i + my_group * cores_per_group;
  }
  for (int i = 0; i < n_group; i++) {
    root_ranks[i] = i * cores_per_group;
  }

  MPI_Comm_group(MPI_COMM_WORLD, &gl_group);
  MPI_Group_incl(gl_group, cores_per_group, ranks, &group);
  MPI_Group_incl(gl_group, n_group, root_ranks, &root_group);
  MPI_Comm_create(MPI_COMM_WORLD, group, &group_comm);
  MPI_Comm_create(MPI_COMM_WORLD, root_group, &root_comm);
  int idx_beg = my_group * arg_size;
#else
  MPI_Comm group_comm = 0;
  MPI_Comm root_comm = 0;
  int my_group = 0;
#endif

  double Lx = atof(argv[1]);
  double Ly = Lx;
  double D_0 = atof(argv[2]);
  double omega0 = atof(argv[3]);
  double rho_0 = atof(argv[4]);
  unsigned long long seed = atoi(argv[5]);
  int n_step = atof(argv[6]);
  std::string ini_mode = argv[7];

  Vec_2<double> gl_l(Lx, Ly);;
  int gl_par_num = int(gl_l.x * gl_l.y * rho_0);
  run(gl_par_num, gl_l, D_0, omega0, seed, n_step, ini_mode, group_comm, root_comm);

#ifdef USE_MPI
  delete[] ranks;
  delete[] root_ranks;
  MPI_Finalize();
#endif
}

#endif
