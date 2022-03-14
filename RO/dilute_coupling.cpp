#include "dilute_coupling.h"
#ifdef DILUTE_COUPLING
#include "run2D.h"
#include "disorder2D.h"

void run_dilute_coupling(int gl_par_num, const Vec_2<double>& gl_l,
                         double eta, double rho_s, double eps,
                         unsigned long long seed, unsigned long long seed2, int n_step,
                         std::string& ini_mode, MPI_Comm group_comm, MPI_Comm root_comm) {
  typedef BiNode<Bird_2> node_t;
  int my_rank = 0;
  int tot_proc = 1;
#ifdef USE_MPI
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
#endif
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;
  const double v0 = 0.5;
  double rho_0 = gl_par_num / (gl_l.x * gl_l.y);

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
    snprintf(basename, 100, "%g_%.3f_%.5f_%.3f_%llu_%03llu", gl_l.x, eta, rho_s, eps, seed, seed2);
  } else {
    snprintf(basename, 100, "%g_%g_%.3f_%.5f_%.3f_%llu_%03llu", gl_l.x, gl_l.y, eta, rho_s, eps, seed, seed2);
  }
  snprintf(scatter_f, 200, "data/scatter_%s.bin", basename);


  // initialize random obstacles
  int n_ob = round(rho_s * gl_l.x * gl_l.y);
  Ranq1 myran(seed + my_rank);
  DiluteScatter disorder(dm, grid, n_ob, myran, group_comm, scatter_f);

  // initialize the location and orientation of particles
  int t_beg;
  ini_particles(gl_par_num, p_arr, ini_mode, myran, seed2, cl, dm, t_beg);

  int ave_t_beg = 400000;
  int ave_t_beg1 = t_beg >= ave_t_beg ? 0 : ave_t_beg - t_beg;
  int ave_dt1 = 100000;
  int field_t_beg1 = 0;
  int field_dt1 = 10000;
  snprintf(logfile, 200, "data/DC_%s_%d.log", basename, t_beg);
  snprintf(phifile, 200, "data/DC_%s_%d.dat", basename, t_beg);
  snprintf(snapfile, 200, "data/snap/DC_%s", basename);
  snprintf(field_f1, 200, "data/DC_field_%s_%d_%d.bin", basename, field_dt1, t_beg);
  snprintf(time_ave_f1, 200, "data/DC_ave_%s_%d_%d.bin", basename, ave_dt1, t_beg);

  exporter::LogExporter log(logfile, 0, n_step * 2, 10000, gl_par_num, group_comm);
  exporter::OrderParaExporter_2 order_ex(phifile, 0, n_step * 2, 100, gl_l, group_comm);
  exporter::SnapExporter snap_ex(snapfile, 0, n_step * 2, 10000, t_beg, group_comm);
  exporter::FeildExporter field_ex1(field_f1, field_t_beg1, n_step * 2, field_dt1, grid, dm, 4);
  exporter::TimeAveFeildExporter ave_ex1(time_ave_f1, ave_t_beg1, n_step * 2, ave_dt1, grid, dm, 4);

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
  double eps2PI = eps * 2.0 * PI;
  Ranq2 myran2(myran.int64() + seed2 + t_beg);
  auto single_move = [eta2PI, v0, &myran2, &dm, &disorder, eps2PI](node_t& p) {
    double noise_strength;
    if (disorder.within(p)) {
      noise_strength = eps2PI;
    } else {
      noise_strength = eta2PI;
    }
    double noise = (myran2.doub() - 0.5) * noise_strength;
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

#ifndef USE_MPI
  auto update_one_step = [&p_arr, &cl, &for_all_pair_force_fast, &single_move](int i) {
    for_all_pair_force_fast();
    integrate(p_arr, cl, single_move);
    out(i, p_arr);
  };
#else
  auto run_one_step = [&p_arr, &cl, &comm, &for_all_pair_force_fast, &single_move, &out](int i) {
    cal_force(p_arr, cl, comm, for_all_pair_force_fast);
    integrate(p_arr, cl, single_move, comm);
    out(i, p_arr);
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

  int arg_size = 9;
  if ((argc - 1) % arg_size != 0) {
    std::cout << "Error, argc = " << argc << std::endl;
    exit(1);
  }
  int n_group = (argc-1)/arg_size;
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

  double Lx = atof(argv[1+idx_beg]);
  double Ly = atof(argv[2+idx_beg]);
  double eta = atof(argv[3+idx_beg]);
  double eps = atof(argv[4+idx_beg]);
  double rho_s = atof(argv[5+idx_beg]);
  unsigned long long seed1 = atoi(argv[6+idx_beg]);
  unsigned long long seed2 = atoi(argv[7+idx_beg]);
  int n_step = atof(argv[8+idx_beg]);
  std::string ini_mode = argv[9+idx_beg];

  Vec_2<double> gl_l(Lx, Ly);
  double rho0 = 1.;
  int gl_par_num = int(gl_l.x * gl_l.y * rho0);
  run_dilute_coupling(gl_par_num, gl_l, eta, rho_s, eps, seed1, seed2, n_step, ini_mode, group_comm, root_comm);

#ifdef USE_MPI
  delete[] ranks;
  delete[] root_ranks;

  MPI_Finalize();
#endif
}
#endif
