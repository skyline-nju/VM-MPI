#pragma once
#include "domain3D.h"
#include "cellList3D.h"

int ini_my_par_num(int gl_par_num);

template <typename TNode, typename TRan>
void ini_rand(std::vector<TNode>& p_arr, int gl_par_num, TRan& myran,
              CellListNode_3<TNode>& cl, Domain_3 &dm) {
  const int my_par_num = ini_my_par_num(gl_par_num);
#ifdef USE_MPI
  p_arr.reserve(my_par_num * 2);
  dm.set_max_buf_size(gl_par_num, 10.);
#else
  p_arr.reserve(my_par_num);
#endif
  for (int i = 0; i < my_par_num; i++) {
    p_arr.emplace_back(myran, dm.l(), dm.origin());
  }
  cl.create(p_arr);
}

template <typename TNode>
void cal_force(std::vector<TNode> &p_arr,
               CellListNode_3<TNode> &cl,
               const Domain_3 &dm) {
#ifdef USE_MPI
  int n_ghost = 0;
  comm_par_before_interact(dm.neighbor, dm.inner_shell, dm.max_buf_size(),
    p_arr, cl, n_ghost);
#endif

  auto f_ij = [&dm](TNode *pi, TNode *pj) {
    pi->interact(*pj, dm);
  };
  cl.for_each_pair(f_ij);

#ifdef USE_MPI
  clear_ghost_after_interact(cl, dm.outer_shell, p_arr, n_ghost);
#endif
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

#ifdef USE_MPI
  comm_par_after_move(dm.neighbor, dm.outer_shell, dm.max_buf_size(), p_arr, cl);
#endif
}

template <typename TNode, typename TInteract, typename TIntegrate>
void run(std::vector<TNode> &p_arr, TInteract interact,
         TIntegrate integrate, int n_step, int sep) {
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  auto t1 = std::chrono::system_clock::now();
  for (int i = 1; i <= n_step; i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    interact(p_arr);
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

