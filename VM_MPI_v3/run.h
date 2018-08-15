#pragma once

#include "domain3D.h"
#include "cellList3D.h"

template <typename TNode, typename TRan>
void ini_rand(std::vector<TNode>& p_arr, int n_par, TRan& myran,
              CellListNode_3<TNode>& cl, const Domain_3 &dm) {
  p_arr.reserve(n_par);
  for (int i = 0; i < n_par; i++) {
    p_arr.emplace_back(myran, dm.gl_l(), dm.origin());
  }
  cl.create(p_arr);
}

template <typename TNode>
void cal_force(const CellListNode_3<TNode> &cl, const Domain_3 &dm) {
  auto f_ij = [&dm](TNode *pi, TNode *pj) {
    pi->interact(*pj, dm);
  };
  cl.for_each_pair(f_ij);
}

template <typename TNode, typename TRan>
void integrate(std::vector<TNode>& p_arr, TRan& myran,
               CellListNode_3<TNode>& cl, const Domain_3 &dm,
               double eta, double v0) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    (*it).move(eta, v0, myran, dm);
  }
  cl.recreate(p_arr);
}

template <typename TNode, typename TInteract, typename TIntegrate>
void run(std::vector<TNode> &p_arr, TInteract interact,
         TIntegrate integrate, int n_step, int sep) {
  auto t1 = std::chrono::system_clock::now();
  for (int i = 1; i <= n_step; i++) {
    interact();
    integrate(p_arr);

    if (i % sep == 0) {
      double phi;
      Vec_3<double> v_mean{};
      cal_order_para(p_arr, phi, v_mean);
      std::cout << phi << "\t" << v_mean << std::endl;
    }
  }
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
  std::cout << n_step * p_arr.size() / elapsed_time.count() << std::endl;
}

