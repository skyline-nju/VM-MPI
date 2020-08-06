#include "config.h"
#include "mpi.h"
#include "communicator2D.h"
#include "domain2D.h"
#include "cellList2D.h"
#include "rand.h"
#include "run2D.h"
#include "particle2D.h"
#ifdef OUTPUT_ON
#include "exporter2D.h"
#endif

using namespace std;
#ifndef REPLICAS
int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  double Lx = atof(argv[1]);
  double Ly = atof(argv[2]);
  double eta = atof(argv[3]);
  double eps = atof(argv[4]);
  int seed = atoi(argv[5]);
  int n_step = atoi(argv[6]);
  Vec_2<double> gl_l(Lx, Ly);

  double rho0 = 1;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  //Ranq1 myran(seed + my_rank);
  double v0 = 0.5;
  double r_cut = 1.;
  int gl_par_num = int(gl_l.x *gl_l.y * rho0);

  double theta0 = 0.;
  if (argc == 7) {
    run_quenched_ini_ordered(gl_par_num, gl_l, eta, eps, seed, n_step, theta0, MPI_COMM_WORLD);
  } else if (argc == 8) {
    int cores_per_group = atoi(argv[7]);
    int tot_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
    int n_group = tot_proc / cores_per_group;
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
    unsigned long long my_seed = seed + my_group;
    run_quenched_ini_ordered(gl_par_num, gl_l, eta, eps, my_seed, n_step, theta0, group_comm, root_comm);
    delete[] ranks;
    delete[] root_ranks;
  }

  MPI_Finalize();
}

#else

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  double Lx = 256;
  double Ly = 256;
  double eta = 0.18;
  int seed = atoi(argv[1]);
  int n_step = 100000;
  Vec_2<double> gl_l(Lx, Ly);

  double rho0 = 1.;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  double v0 = 0.5;
  double r_cut = 1.;
  int gl_par_num = int(gl_l.x * gl_l.y * rho0);

  int cores_per_group = 4;
  int tot_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  int n_group = tot_proc / cores_per_group;
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

  double eps = atof(argv[my_group * 2 + 2]);
  double theta0 = atof(argv[my_group * 2 + 3]);
  run_quenched_ini_ordered(gl_par_num, gl_l, eta, eps, seed, n_step, theta0, group_comm, root_comm);
  delete[] ranks;
  delete[] root_ranks;
  MPI_Finalize();
}
#endif
