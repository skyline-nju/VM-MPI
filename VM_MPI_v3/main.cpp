#include "singleDomain3D.h"

using namespace std;

int main(int argc, char* argv[]) {
  double L = 48;
  Vec_3<double> l(L, L, L);
  double eta = 0.2;
  double rho0 = 1.;
  unsigned long long seed = 12;
  int n_step = 10000;
  double v0 = 0.5;

  SingleDomain_3<UniNode<VicsekPar_3>> domain(l, rho0, seed);
  domain.ini_rand();
  domain.eval_elapsed_time(eta, v0, n_step, 200);
};