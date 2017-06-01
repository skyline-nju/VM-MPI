#include <iostream>
#include <vector>
#include "mpi.h"
#include "domain.h"

using namespace std;

int main(int argc, char *argv[]) {
  double Lx = 180;
  double Ly = atof(argv[1]);
  //double Ly = 600;
  double eta = 0.35;
  double rho0 = 1;

  MPI_Init(&argc, &argv);
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  double t_beg;

  StaticDomain subdomain(Lx, Ly, size, rank, 123, eta, 0, rho0);
  subdomain.create_from_snap("s_0.35_0_1_180_100_123_01000000.bin");
  //subdomain.create_particle_random(int(rho0 * Lx * Ly / size));
  if (rank == 0) {
    t_beg = MPI_Wtime();
  }
  for (int i = 0; i < 1000; i++) {
    subdomain.one_step_MPI(eta, i);
    if (i % 1000 == 0) {
      cout << "t = " << i << endl;
    }
  }

  if (rank == 0) {
    cout << "elapsed time = " << MPI_Wtime() - t_beg << endl;
    //Cell a[10];
    //cout << sizeof(&a[0]) << endl;
    //Cell *p = &a[9];
    //int x = p - &a[0];
    //cout << x << endl;
  }
  MPI_Finalize();
  return 0;
}