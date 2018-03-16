#ifndef SUBDOMAIN_H
#define SUBDOMAIN_H
#include <vector>
#include "particle.h"
#include "vect.h"
#include "cellList.h"
#include "boundary.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

typedef unsigned long long int Ull;

class DomainBase {
public:
  DomainBase(double Lx, double Ly, double Rho, double eta0,
    Ull seed0, double v = 0.5);

  virtual ~DomainBase();

  virtual void ini_rand() = 0;

  virtual void cal_force() = 0;

  virtual void integrate() = 0;

  virtual void integrate2() = 0;

  virtual void cal_order_para(double &phi, double &theta) const = 0;

protected:
  Ran *myran;
  Vec_2<double> system_size;
  double rho0;
  double eta;
  double v0;
  Ull seed;
  int nPar;
};

#ifndef USE_MPI

class SingleDomain : public DomainBase {
public:
  SingleDomain(double Lx, double Ly, double Rho, double eta0, Ull seed0);

  ~SingleDomain() { delete pbc; delete cl; }

  void ini_rand();

  void cal_force();

  void integrate();

  void integrate2() {}

  void cal_order_para(double &phi, double &theta) const;

protected:
  PBC_2 *pbc;
  CellListIdx_2 *cl;
  std::vector<Par1> p_arr;
};

class SingleDomain2 : public DomainBase {
public:
  SingleDomain2(double Lx, double Ly, double Rho, double eta0, Ull seed0);

  ~SingleDomain2() { delete pbc; delete cl; }

  void ini_rand();

  void cal_force();

  void integrate();

  void integrate2();

  void cal_order_para(double &phi, double &theta) const;

protected:
  PBC_2 * pbc;
  CellListNode_2 *cl;
  std::vector<Node> p_arr;

};

class SingleDomain3 : public DomainBase {
public:
  SingleDomain3(double Lx, double Ly, double Rho, double eta0, Ull seed0);

  ~SingleDomain3() { delete pbc; delete cl; }

  void ini_rand();

  void cal_force();

  void integrate();

  void integrate2();

  void cal_order_para(double &phi, double &theta) const;

protected:
  PBC_2 * pbc;
  CellListBiNode_2 *cl;
  std::vector<BiNode> p_arr;

};

#else
class StaticSubDomain : public DomainBase {
public:
  StaticSubDomain(double Lx, double Ly, double Rho, double eta0, Ull seed0);

  ~StaticSubDomain() { delete pbc_x; delete cl;}

  void ini_rand();

  void cal_force();

  void integrate();


protected:
  CellListIdx_2 * cl;
  PBC_X_2 * pbc_x;
  int my_rank;
  int pre_rank;
  int next_rank;
  int tot_processor;
  int my_nPar;
};
#endif

#endif


