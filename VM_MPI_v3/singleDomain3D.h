#pragma once

#include "cellList3D.h"
#include "rand.h"
#include "boundary3D.h"
#include "particle3D.h"

template <typename TNode>
class SingleDomain_3 {
public:
  typedef unsigned long long int Ull;
  typedef Vec_3<double> Vec3d;
  SingleDomain_3(const Vec3d &l, double rho0, Ull seed, double r_cut = 1.);

  void ini_rand();
  void cal_force();
  void integrate(double eta, double v0 = 0.5);

  void eval_elapsed_time(double eta, double v0, int n_step, int sep);

private:
  Ranq1 myran_;
  Vec_3<double> l_;
  CellListNode_3<TNode> cell_list_;
  PBC_xyz_3 bc_;
  std::vector<TNode> p_arr_;

  int n_par_;
};

template<typename TNode>
SingleDomain_3<TNode>::SingleDomain_3(const Vec3d & l, double rho0,
                                      Ull seed, double r_cut)
  : myran_(seed), l_(l), cell_list_(l, 1.), bc_(l){
  n_par_ = int(l.x * l.y * l.z * rho0);
  p_arr_.reserve(n_par_);
}

template <typename TNode>
void SingleDomain_3<TNode>::ini_rand() {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran_, l_, Vec3d());
  }
  cell_list_.create(p_arr_);
  double phi;
  Vec3d v_mean;
  cal_order_para(p_arr_, phi, v_mean);
  std::cout << "t = 0\t phi = " << phi << "\t" << v_mean << std::endl;

}

template <typename TNode>
void SingleDomain_3<TNode>::cal_force() {
  auto f_ij = [this](TNode *pi, TNode *pj) {
    pi->interact(*pj, bc_);
  };
  cell_list_.for_each_pair(f_ij);
}

template <typename TNode>
void SingleDomain_3<TNode>::integrate(double eta, double v0) {
  const auto end = p_arr_.end();
  for (auto it = p_arr_.begin(); it != end; ++it) {
    (*it).move(eta, v0, myran_, bc_);
  }
  cell_list_.recreate(p_arr_);
}

template <typename TNode>
void SingleDomain_3<TNode>::eval_elapsed_time(double eta, double v0,
                                              int n_step, int sep) {
  auto t1 = std::chrono::system_clock::now();
  for (int i = 1; i <= n_step; i++) {
    cal_force();
    integrate(eta, v0);

    if (i % sep == 0) {
      double phi;
      Vec3d v_mean;
      cal_order_para(p_arr_, phi, v_mean);
      std::cout << phi << "\t" << v_mean << std::endl;
    }
  }
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
  std::cout << n_step * p_arr_.size() / elapsed_time.count() << std::endl;
}
