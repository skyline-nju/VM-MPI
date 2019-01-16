#pragma once
#include "config.h"
#include "domain2D.h"
#include "cellList2D.h"
#include "disorder2D.h"
#include "particle2D.h"
#include "exporter2D.h"
#include <iomanip>
#include "communicator2D.h"

int ini_my_par_num(int gl_par_num);

template <typename TNode, typename TRan>
void ini_rand(std::vector<TNode>& p_arr, int gl_par_num, TRan& myran,
              CellListNode_2<TNode>& cl, Domain_2 &dm) {
  const int my_par_num = ini_my_par_num(gl_par_num);
#ifdef USE_MPI
  p_arr.reserve(my_par_num * 3);
#else
  p_arr.reserve(my_par_num);
#endif
  for (int i = 0; i < my_par_num; i++) {
    p_arr.emplace_back(myran, dm.l(), dm.origin());
  }
  cl.create(p_arr);
}

template <typename TNode>
void cal_force(std::vector<TNode> &p_arr, CellListNode_2<TNode> &cl, Communicator &comm) {
  int n_ghost = 0;
  comm.comm_before_cal_force(p_arr, cl, n_ghost);

  auto f1 = [](TNode *pi, TNode *pj) {
    pi->interact(*pj);
  };
  auto f3 = [](TNode *pi, TNode *pj, const Vec_2<double> &offset) {
    pi->interact(*pj, offset);
  };
  cl.for_each_pair_fast(f1, f3);

  comm.clear_padded_particles(cl, p_arr, n_ghost);
}

template <typename TNode, typename UniFunc, typename T1, typename T2>
void integrate(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, UniFunc f1,
               Communicator &comm, std::vector<T1> &n_arr, std::vector<Vec_2<T2>> &v_arr) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    f1(*it);
  }
  cl.recreate(p_arr, n_arr, v_arr);
  comm.comm_after_integration(p_arr, cl, n_arr, v_arr);
}

template <typename TNode, typename TInteract, typename TIntegrate, typename TOut>
void run(std::vector<TNode> &p_arr, TInteract interact_all,
         TIntegrate integrate_all, TOut out, int n_step) {
  for (int i = 1; i <= n_step; i++) {
    interact_all(p_arr);

    integrate_all(p_arr);

    out(i, p_arr);
  }
}

void run(int gl_par_num, const Vec_2<double> &gl_l,
         double eta, double eps, unsigned long long seed,
         int n_step, int n_save_snap, int block_size);

void run_birth_death(int gl_par_num, const Vec_2<double>& gl_l,
                     double eta, double alpha, unsigned long long seed,
                     int n_step, int snap_interval, int snap_block_size,
                     int box_len_birth);
