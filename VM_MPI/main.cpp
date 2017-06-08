#include <iostream>
#include <vector>
#include "mpi.h"
#include "domain.h"
#include "cmdline.h"

using namespace std;

int main(int argc, char *argv[]) {
  cmdline::parser cmd;
  cmd.add<double>("eta", '\0', "noise strength", true);
  cmd.add<double>("eps", '\0', "disorder strength", false, 0);
  cmd.add<double>("rho0", '\0', "density of particles", true);
  cmd.add<double>("Lx", '\0', "Lx", true);
  cmd.add<double>("Ly", '\0', "Ly", true);
  cmd.add<int>("nstep", 't', "total steps to run", true);
  cmd.add<unsigned long long int>("seed", 's', "seed of random number", true);
  cmd.add<string>("fin", 'f', "input snapshot", false, "");
  cmd.add("dynamic", '\0', "dynamic domain decomposition");
  cmd.add<int>("rate", '\0', "refresh rate of domain", false, 1);
  cmd.add<int>("cg_dt", '\0', "interval of coarse-grained snapshots", false, 0);
  cmd.add<double>("cg_exp", '\0',
                  "exponent to output coarse-grained snapshots in log scales",
                  false, 0.);
  cmd.parse_check(argc, argv);
  double eta = cmd.get<double>("eta");
  double eps = cmd.get<double>("eps");
  double rho0 = cmd.get<double>("rho0");
  double Lx = cmd.get<double>("Lx");
  double Ly = cmd.get<double>("Ly");
  unsigned long long seed = cmd.get<unsigned long long>("seed");
  int nStep = cmd.get<int>("nstep");
  int rate = cmd.get<int>("rate");

  MPI_Init(&argc, &argv);
  BasicDomain *subdomain;
  double magnification;
  if (cmd.exist("dynamic")) {
    subdomain = new DynamicDomain(cmd);
    magnification = 1.5;
  }
  else {
    subdomain = new StaticDomain(cmd);
    magnification = 3.0;
  }
  if (cmd.exist("fin"))
    subdomain->create_from_snap(cmd, magnification);
  else
    subdomain->create_particle_random(int(Lx * Ly * rho0), magnification);

  for (int i = 1; i <= nStep; i++) {
    subdomain->one_step(i);
  }
  delete subdomain;
  MPI_Finalize();
  return 0;
}