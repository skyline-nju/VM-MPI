#pragma once

#include "vect.h"
#include "comn.h"
#include <vector>

class Domain_3 {
public:
  typedef Vec_3<double> Vec3d;

  explicit Domain_3(const Vec3d &l, const Vec3d &origin = Vec3d());

  template <typename TNode, typename TRan, typename TCellList>
  void ini_rand(std::vector<TNode> &p_arr, int n_par, TRan &myran, TCellList &cl);

  template <typename TNode, typename TCellList>
  void cal_pair_force(std::vector<TNode> &p_arr, TCellList &cl);

  template <typename TNode, typename TRan, typename TCellList>
  void integrate(std::vector<TNode> &p_arr, TRan &myran, TCellList &cl,
                 double eta, double v0 = 0.5);

  template <typename TNode, typename TRan, typename TCellList>
  void eval_elapsed_time(std::vector<TNode> &p_arr, TRan &myran, TCellList &cl,
                         double eta, double v0, int n_step, int sep);
protected:
  Vec3d l_;
  Vec3d half_l_;
  Vec3d origin_;
  Vec3d end_pnt_;
};

template <typename TNode, typename TRan, typename TCellList>
void Domain_3::ini_rand(std::vector<TNode> &p_arr, int n_par,
                        TRan &myran, TCellList &cl) {
  p_arr.reserve(n_par);
  for (int i = 0; i < n_par; i++) {
    p_arr.emplace_back(myran, l_, Vec3d());
  }
  cl.create(p_arr);
}

template <typename TNode, typename TCellList>
void Domain_3::cal_pair_force(std::vector<TNode> &p_arr, TCellList &cl) {
  auto f_ij = [this](TNode *pi, TNode *pj) {
    pi->interact(*pj, l_, half_l_);
  };
  cl.for_each_pair(f_ij);
}

template <typename TNode, typename TRan, typename TCellList>
void Domain_3::integrate(std::vector<TNode>& p_arr, TRan& myran, TCellList &cl,
                         double eta, double v0) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    (*it).move(eta, v0, myran, origin_, end_pnt_, l_);
  }
  cl.recreate(p_arr);
}

template <typename TNode, typename TRan, typename TCellList>
void Domain_3::eval_elapsed_time(std::vector<TNode>& p_arr, TRan& myran, TCellList &cl,
                                 double eta, double v0, int n_step, int sep) {
  auto t1 = std::chrono::system_clock::now();
  for (int i = 1; i <= n_step; i++) {
    cal_pair_force(p_arr, cl);
    integrate(p_arr, myran, cl, eta, v0);

    if (i % sep == 0) {
      double phi;
      Vec3d v_mean;
      cal_order_para(p_arr, phi, v_mean);
      std::cout << phi << "\t" << v_mean << std::endl;
    }
  }
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
  std::cout << n_step * p_arr.size() / elapsed_time.count() << std::endl;
}


class ExtDomain_3: public Domain_3 {
public:
  ExtDomain_3(const Vec3d &l, const Vec3d &origin,
              const Vec_3<bool> &flag_ext);

protected:
  Vec_3<bool> flag_ext_;
};


