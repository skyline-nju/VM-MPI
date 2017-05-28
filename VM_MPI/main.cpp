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

  SubDomain subdomain(Lx, Ly, size, rank, 123);
  subdomain.create_particle_random(int(rho0 * Lx * Ly / size));
  for (int i = 0; i < 10000; i++) {
    subdomain.one_step_MPI(eta, i);
  }
  MPI_Finalize();
}