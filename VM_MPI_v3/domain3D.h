#pragma once

#include <vector>
#include "vect.h"
#include "comn.h"
#define USE_MPI

#ifdef USE_MPI
#include "communicator3D.h"
#endif

class Domain_3 {
public:
  typedef Vec_3<double> Vec3d;

  Domain_3() = default;

  explicit Domain_3(const Vec3d &l, const Vec3d &origin = Vec3d())
    : l_(l), half_l_(0.5 * l), origin_(origin), end_pnt_(origin + l) {}

  const Vec3d & l() const { return l_; }
  const Vec3d & half_l() const { return half_l_; }
  const Vec3d & origin() const { return origin_; }
  const Vec3d & end_pnt() const { return end_pnt_; }


  template <typename TNode, typename TRan, typename TCellList>
  void ini_rand(std::vector<TNode> &p_arr, int n_par, TRan &myran, TCellList& cl);

  template <typename TNode, typename TCellList>
  void cal_pair_force(std::vector<TNode> &p_arr,
                      const TCellList &cl) const;

  template <typename TNode, typename TRan, typename TCellList>
  void integrate(std::vector<TNode> &p_arr, TRan &myran, TCellList &cl,
                 double eta, double v0 = 0.5);

  static Vec_3<int> partitioning();
protected:
  Vec3d l_;
  Vec3d half_l_;
  Vec3d origin_;
  Vec3d end_pnt_;
};

template <typename TNode, typename TRan, typename TCellList>
void Domain_3::ini_rand(std::vector<TNode>& p_arr, int n_par, TRan& myran, TCellList& cl) {
  p_arr.reserve(n_par);
  for (int i = 0; i < n_par; i++) {
    p_arr.emplace_back(myran, l_, Vec3d());
  }
  cl.create(p_arr);
}

template<typename TNode, typename TCellList>
void Domain_3::cal_pair_force(std::vector<TNode>& p_arr,
                             const TCellList & cl) const {
  auto f_ij = [this](TNode *pi, TNode *pj) {
    pi->interact(*pj, *this);
  };
  cl.for_each_pair(f_ij);
}

template <typename TNode, typename TRan, typename TCellList>
void Domain_3::integrate(std::vector<TNode>& p_arr, TRan& myran, TCellList& cl,
                        double eta, double v0) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    (*it).move(eta, v0, myran, *this);
  }
  cl.recreate(p_arr);
}

template <typename TNode, typename TInteract, typename TIntegrate>
void run(std::vector<TNode> &p_arr, TInteract interact,
         TIntegrate integrate, int n_step, int sep) {
  auto t1 = std::chrono::system_clock::now();
  for (int i = 1; i <= n_step; i++) {
    interact(p_arr);
    integrate(p_arr);

    if (i % sep == 0) {
      double phi;
      Vec_3<double> v_mean;
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
  ExtDomain_3(const Vec_3<int> &domain_size, const Vec3d &gl_l, double r_cut, int gl_np);

  const Vec_3<int> & cell_size() const { return my_cell_size_; }
  const Vec_3<double> & cell_len() const { return cell_len_; }
  const Vec_3<bool> &flag_ext() const { return flag_ext_; }

  template <typename TNode, typename TRan, typename TCellList>
  void ini_rand(std::vector<TNode> &p_arr, TRan &myran, TCellList& cl);

  template <typename TNode, typename TCellList>
  void cal_pair_force(std::vector<TNode> &p_arr,
    const TCellList &cl) const;

  template <typename TNode, typename TRan, typename TCellList>
  void integrate(std::vector<TNode> &p_arr, TRan &myran, TCellList &cl,
    double eta, double v0 = 0.5);
protected:
  Vec_3<bool> flag_ext_;
  Vec_3<int> size_;
  Vec_3<int> rank_;
  Vec_3<int> gl_cell_size_;
  Vec_3<int> my_cell_size_;
  Vec_3<int> my_first_cell_;
  Vec_3<double> cell_len_;


  int my_par_num_;
  int gl_par_num_;
  int my_proc_;
  int proc_size_;
  int neighbor_proc_[3][2]{};
};

template <typename TNode, typename TRan, typename TCellList>
void ExtDomain_3::ini_rand(std::vector<TNode>& p_arr, TRan& myran, TCellList& cl) {
  if (proc_size_ == 1)
    p_arr.reserve(gl_par_num_);
  else
    p_arr.reserve(gl_par_num_ / proc_size_ * 2);
  for (int i = 0; i < my_par_num_; i++) {
    p_arr.emplace_back(myran, l_, origin_);
  }
  cl.create(p_arr);
}

template <typename TNode, typename TCellList>
void ExtDomain_3::cal_pair_force(std::vector<TNode>& p_arr,
                                 const TCellList& cl) const {
  auto f_ij = [this](TNode *pi, TNode *pj) {
    pi->interact(*pj, *this);
  };
  cl.for_each_pair(f_ij);
}

template <typename TNode, typename TRan, typename TCellList>
void ExtDomain_3::integrate(std::vector<TNode>& p_arr, TRan& myran, TCellList& cl,
                            double eta, double v0) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    (*it).move(eta, v0, myran, *this);
  }
  cl.recreate(p_arr);
}



