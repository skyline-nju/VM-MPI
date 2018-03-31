#include <iostream>
#include "subDomain.h"
#include <ctime>
#include <chrono>
#include "singleDomain2D.h"
#include "node.h"

using namespace std;

//template <typename T>
//void foo(T x, vector<T*> &y) {
//	cout << "a" << endl;
//}
//
//template <typename T>
//void foo(T x, vector<T> &y) {
//	cout << "b" << endl;
//}

int main(int argc, char* argv[]) {

	//std::vector<double *> vect_a;
	//std::vector<double> vect_b;

	//double scalar = 1;
	//foo(scalar, vect_a);
	//foo(scalar, vect_b);




#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  //int my_rank;
  //int size;
  //MPI_Comm_size(MPI_COMM_WORLD, &size);
  //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  //cout << my_rank << " of " << size << " procs\n";

  double L = 64;
  double rho0 = 1;
  double eta = 0.1;
  unsigned long long seed = 1;
  int n_step = 10000;
  int t_start = 4000;
  int t_sep = 100;
#ifndef USE_MPI
  SingleDomain2 domain(L, L, rho0, eta, seed);
  domain.ini_rand();
  double phi_mean = 0;
  int count = 0;
  double phi, theta;
  auto t1 = std::chrono::system_clock::now();
  for (int i = 0; i < n_step; i++) {
    domain.cal_force();
    domain.integrate();
    if (i >= t_start && i % t_sep == 0) {
      domain.cal_order_para(phi, theta);
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

  //delete domain;
  //domain = new SingleDomain(L, L, rho0, eta, seed);
  //domain->ini_rand();

  //phi_mean = 0;
  //count = 0;
  //t1 = std::chrono::system_clock::now();
  //for (int i = 0; i < n_step; i++) {
  //  domain->cal_force();
  //  if (i >= t_start && i % t_sep == 0) {
  //    double phi, theta;
  //    domain->cal_order_para(phi, theta);
  //    //std::cout << phi << "\t" << theta << "\n";
  //    phi_mean += phi;
  //    count++;
  //  }
  //  domain->integrate();
  //}
  //t2 = std::chrono::system_clock::now();
  //elapsed_time = t2 - t1;
  //std::cout << elapsed_time.count() << std::endl;
  //std::cout << "phi = " << phi_mean / count << std::endl << std::endl;

  //SDomain_2<Par1> sd(L, L, rho0, seed);
  //sd.ini_rand();
  //phi_mean = 0;
  //count = 0;
  //t1 = std::chrono::system_clock::now();
  //for (int i = 0; i < n_step; i++) {
  //  sd.cal_force();
  //  if (i >= t_start && i % t_sep == 0) {
  //    double phi, theta;
  //    sd.cal_order_para(phi, theta);
  //    //std::cout << phi << "\t" << theta << "\n";
  //    phi_mean += phi;
  //    count++;
  //  }
  //  sd.integrate(eta);
  //}
  //t2 = std::chrono::system_clock::now();
  //elapsed_time = t2 - t1;
  //std::cout << elapsed_time.count() << std::endl;
  //std::cout << "phi = " << phi_mean / count << std::endl;

  //SDomain_2<Par1, std::list<Par1*>> sd2(L, L, rho0, seed);
  SDomain_2<UniNode<Par1>, UniNode<Par1>*> sd2(L, L, rho0, seed);

  sd2.ini_rand();
  phi_mean = 0;
  count = 0;
  t1 = std::chrono::system_clock::now();
  for (int i = 0; i < n_step; i++) {
    sd2.cal_force();
    sd2.integrate(eta);
    if (i >= t_start && i % t_sep == 0) {
      double phi, theta;
      sd2.cal_order_para(phi, theta);
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