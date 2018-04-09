#include <iostream>
#include <ctime>
#include <chrono>
//#include "singleDomain2D.h"
#include "node.h"
#include "particle.h"
#include "singleDomain2D.h"
//#define USE_MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

int main(int argc, char* argv[]) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  int my_rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

#ifndef USE_MPI
  double L = 64;
  double rho0 = 1;
  double eta = 0.1;
  unsigned long long seed = 1;
  int n_step = 10000;
  int t_start = 2000;
  int t_sep = 100;

  {
    Single_domain_2<BiNode<Par1>> s(L, L, rho0, seed);
    s.ini_rand();
    s.eval_elapsed_time(eta, n_step, t_start, t_sep, 1);
  }

  {
    Single_domain_2<BiNode<Par1>> s(L, L, rho0, seed);
    s.ini_rand();
    s.eval_elapsed_time(eta, n_step, t_start, t_sep, 3);
  }

  {
    Single_domain_2<BiNode<Par1>> s(L, L, rho0, seed);
    s.ini_rand();
    s.eval_elapsed_time(eta, n_step, t_start, t_sep, 2);
  }

#else
  MPI_Finalize();
#endif

}