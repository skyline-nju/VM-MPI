#ifndef SINGLE_DOMAIN_2_H
#define SIGNLE_DOMAIN_2_H
#include "cellList2D.h"
#include "rand.h"
#include "boundary.h"
#include "particle.h"

template <typename _TNode>
class Single_domain_2 {
public:
  typedef unsigned long long int Ull;
  Single_domain_2(double Lx, double Ly, double Rho, Ull seed0, double r_cut = 1);

  void ini_rand();
  void cal_force();
  void integrate(double eta, double v0 = 0.5);
  void integrate2(double eta, double v0 = 0.5);
  void integrate3(double eta, double v0 = 0.5);
  void cal_order_para(double &phi, double &theta) const;
  void eval_elapsed_time(double eta, int n_steps, int t_start, int t_sep, int mode = 1);
private:
  Ran myran;
  Cell_list_2<_TNode> cell_list;
  PBC_2 pbc;
  std::vector<_TNode> p_arr;
  Vec_2<double> domain_size;

  int nPar;
};

template<typename _TNode>
Single_domain_2<_TNode>::Single_domain_2(double Lx, double Ly, double Rho, Ull seed0, double r_cut):
    myran(seed0), cell_list(Lx, Ly, 0, 0, r_cut), pbc(Lx, Ly) {
  nPar = int(Lx * Ly * Rho);
  p_arr.reserve(nPar);
  domain_size.x = Lx;
  domain_size.y = Ly;
}

template<typename _TNode>
void Single_domain_2<_TNode>::ini_rand() {
  for (int i = 0; i < nPar; i++) {
    p_arr.emplace_back(myran, domain_size.x, domain_size.y);
  }
  cell_list.create(p_arr);
  double phi, theta;
  cal_order_para(phi, theta);
  std::cout << "t = 0\t phi = " << phi << ", theta = " << theta << "\n";
}

template<typename _TNode>
void Single_domain_2<_TNode>::cal_force() {
  auto f_ij = [this](_TNode *pi, _TNode *pj) {
    pi->interact(*pj, pbc);
  };
  cell_list.for_each_pair(f_ij);
}

template<typename _TNode>
void Single_domain_2<_TNode>::integrate(double eta, double v0) {
  for (int i = 0; i < nPar; i++) {
    p_arr[i].move(eta, v0, myran, pbc);
  }
  cell_list.recreate(p_arr);
}

template<typename _TNode>
void Single_domain_2<_TNode>::integrate2(double eta, double v0) {
  auto move = [this, eta, v0](_TNode *p) {
    p->move(eta, v0, myran, pbc);
  };
  cell_list.update(move);
}

template<typename _TNode>
void Single_domain_2<_TNode>::integrate3(double eta, double v0) {
  auto move = [this, eta, v0](_TNode *p) {
    p->move(eta, v0, myran, pbc);
  };
  cell_list.update_by_row(move);
}

template<typename _TNode>
inline void Single_domain_2<_TNode>::cal_order_para(double & phi, double & theta) const {
  func_order_para(p_arr, phi, theta);
}

template<typename _TNode>
inline void Single_domain_2<_TNode>::eval_elapsed_time(double eta, int n_steps, int t_start, int t_sep, int mode) {
  auto t1 = std::chrono::system_clock::now();
  double phi_mean = 0;
  int count = 0;
  for (int i = 0; i < n_steps; i++) {
    cal_force();
    if (mode == 1)
      integrate(eta);
    else if (mode == 3)
      integrate2(eta);
    else if (mode == 2)
      integrate3(eta);
    if (i >= t_start && i % t_sep == 0) {
      double phi, theta;
      cal_order_para(phi, theta);
      phi_mean += phi;
      count++;
    }
  }
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
  std::cout << "phi: " << phi_mean / count << std::endl;
  std::cout << std::endl;
}

#endif