#include <iostream>
#include <vector>
#include "mpi.h"
#include "domain.h"

using namespace std;

int main(int argc, char *argv[]) {
  double Lx = 100;
  double Ly = 300;
  double eta = 0.1;
  double rho0 = 1;

  MPI_Init(&argc, &argv);
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  double t_beg;
  if (rank == 0) {
    t_beg = MPI_Wtime();
  }

  SubDomain subdomain(Lx, Ly, size, rank, 123);
  subdomain.create_particle_random(int(rho0 * Lx * Ly / size));
  double sum_phi = 0;
  int count = 0;
  for (int i = 0; i < 1000000; i++) {
    subdomain.one_step_MPI(eta, i, sum_phi, count);
  }
  if (rank == 0) {
    cout << "elapsed time = " << MPI_Wtime() - t_beg << endl;
  }
  MPI_Finalize();
}