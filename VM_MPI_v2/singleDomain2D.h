#ifndef SINGLEDOMAIN2D_H
#define SINGLEDOMAIN2D_H
#include "cellList2D.h"
#include "rand.h"
#include "boundary.h"

template <typename _TPar, typename _TCell>
class SDomain_2 {
public:
  typedef unsigned long long int Ull;
  SDomain_2(double Lx, double Ly, double Rho, Ull seed0, double r_cut=1);

  void ini_rand();
  void cal_force();
  void integrate(double eta, double v0 = 0.5);
  void integrate2(double eta, double v0 = 0.5);
  void cal_order_para(double &phi, double &theta) const;

private:
  Ran myran;
	Cell_list_2<_TCell> cell_list;
  PBC_2 pbc;
  std::vector<_TPar> p_arr;
  Vec_2<double> domain_size;

  int nPar;
};



template<typename _TPar, typename _TCell>
SDomain_2<_TPar, _TCell>::SDomain_2(double Lx, double Ly, double Rho, Ull seed0, double r_cut):
  myran(seed0), cell_list(Lx, Ly, 0, 0, r_cut), pbc(Lx, Ly) {
  nPar = int(Lx * Ly * Rho);
  p_arr.reserve(nPar);
  domain_size.x = Lx;
  domain_size.y = Ly;

  std::cout << "seed = " << seed0 << "\n";
  std::cout << "nPar = " << nPar << "\n";
}

template<typename _TPar, typename _TCell>
void SDomain_2<_TPar, _TCell>::ini_rand() {
  for (int i = 0; i < nPar; i++) {
    p_arr.emplace_back(myran, domain_size.x, domain_size.y);
  }
  cell_list.create(p_arr);
  double phi, theta;
  cal_order_para(phi, theta);
  std::cout << "t = 0\t phi = " << phi << ", theta = " << theta << "\n";
}

template<typename _TPar, typename _TCell>
void SDomain_2<_TPar, _TCell>::cal_force() {
  auto f_ij = [this](_TPar *pi, _TPar *pj) {
    pi->interact(*pj, pbc);
  };
  cell_list.for_each_pair(f_ij);
  //int count = 0;
  //auto f_ij = [this, &count](_TPar *pi, _TPar *pj) {
  //  pi->interact(*pj, pbc, count);
  //};
  //cell_list.for_each_pair(f_ij);
  //std::cout << "count = " << count << std::endl;
}

template<typename _TPar, typename _TCell>
void SDomain_2<_TPar, _TCell>::integrate(double eta, double v0) {
  for (int i = 0; i < nPar; i++) {
    p_arr[i].move(eta, v0, myran, pbc);
  }
  cell_list.recreate(p_arr);
}

template<typename _TPar, typename _TCell>
void SDomain_2<_TPar, _TCell>::integrate2(double eta, double v0) {
  auto move = [this, eta, v0](_TPar *p) {
    p->move(eta, v0, myran, pbc);
  };
  cell_list.update(move);
}

template<class _TPar, typename _TCell>
inline void SDomain_2<_TPar, _TCell>::cal_order_para(double & phi, double & theta) const {
  func_order_para(p_arr, phi, theta);
}

#endif
