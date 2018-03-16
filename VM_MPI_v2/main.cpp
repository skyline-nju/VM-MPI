#include <iostream>
#include "subDomain.h"
#include <ctime>
#include <chrono>

using namespace std;

int main(int argc, char* argv[]) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  //int my_rank;
  //int size;
  //MPI_Comm_size(MPI_COMM_WORLD, &size);
  //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  //cout << my_rank << " of " << size << " procs\n";

  double L = 100;
  double rho0 = 1;
  double eta = 0.1;
  unsigned long long seed = 1;
  DomainBase *domain;
  int n_step = 20000;
  int t_start = 10000;
  int t_sep = 100;
#ifndef USE_MPI
  domain = new SingleDomain2(L, L, rho0, eta, seed);
  domain->ini_rand();

  auto t1 = std::chrono::system_clock::now();
  double phi_mean = 0;
  int count = 0;
  for (int i = 0; i < n_step; i++) {
    domain->cal_force();
    domain->integrate();
    if (i > t_start && i % t_sep == 0) {
      double phi, theta;
      domain->cal_order_para(phi, theta);
      //std::cout << phi << "\t" << theta << "\n";
      phi_mean += phi;
      count++;
    }
  }
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << elapsed_time.count() << std::endl;
  std::cout << "phi = " << phi_mean / count << std::endl;

  std::cout << std::endl;

  delete domain;
  domain = new SingleDomain3(L, L, rho0, eta, seed);
  domain->ini_rand();

  phi_mean = 0;
  count = 0;
  t1 = std::chrono::system_clock::now();
  for (int i = 0; i < n_step; i++) {
    domain->cal_force();
    domain->integrate2();
    if (i > t_start && i % t_sep == 0) {
      double phi, theta;
      domain->cal_order_para(phi, theta);
      //std::cout << phi << "\t" << theta << "\n";
      phi_mean += phi;
      count++;
    }
  }
  t2 = std::chrono::system_clock::now();
  elapsed_time = t2 - t1;
  std::cout << elapsed_time.count() << std::endl;
  std::cout << "phi = " << phi_mean / count << std::endl;


#else

  MPI_Finalize();
#endif
}