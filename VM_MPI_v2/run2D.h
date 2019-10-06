#pragma once
#include "config.h"
#include "domain2D.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "exporter2D.h"
#include <iomanip>
#include "communicator2D.h"
#ifdef RANDOM_TORQUE
#include "disorder2D.h"
#endif

int ini_my_par_num(int gl_par_num);

void set_particle_num(int gl_par_num, int& my_par_num, int& my_par_num_max);

template <typename TNode, typename TRan>
void ini_rand(std::vector<TNode>& p_arr, int gl_par_num, TRan& myran,
              CellListNode_2<TNode>& cl, Domain_2 &dm) {
  int my_par_num, n_max;
  set_particle_num(gl_par_num, my_par_num, n_max);
  p_arr.reserve(n_max);
  for (int i = 0; i < my_par_num; i++) {
    p_arr.emplace_back(myran, dm.l(), dm.origin());
  }
  cl.create(p_arr);
}

template <typename TNode>
void ini_from_snap(std::vector<TNode>& p_arr, int gl_par_num,
                   CellListNode_2<TNode>& cl,
                   Domain_2 &dm, const char* filename) {
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int my_par_num, n_max;
  set_particle_num(gl_par_num, my_par_num, n_max);

  std::string basename = split(filename, "/").back();
  std::vector<std::string> str_arr = split(basename, "_");
  double lx_in, ly_in, rho0;
  str_to_num(str_arr[1], lx_in);
  str_to_num(str_arr[2], ly_in);
  str_to_num(str_arr[4], rho0);
  int n_in = int(lx_in * ly_in * rho0);

  int n_cols = int(dm.l().x / lx_in);
  int n_rows = int(dm.l().y / ly_in);

  double buf[4];
  std::ifstream fin(filename, std::ios::binary);
  for(int i = 0; i < n_in; i++) {
    fin.read((char*)buf, sizeof(double) * 4);
    for (int row = 0; row < n_rows; row++) {
      for (int col = 0; col < n_cols; col++) {
        p_arr.emplace_back(buf);
        p_arr.back().pos.x += dm.origin().x + col * lx_in;
        p_arr.back().pos.y += dm.origin().y + row * ly_in;
      }
    }
  }

  if (my_rank == 0) {
    std::cout << "ini from " << filename << std::endl;
    std::cout << "n_cols = " << n_cols << std::endl;
    std::cout << "n_rows = " << n_rows << std::endl;

  }
  fin.close();
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

template <typename TNode>
void cal_force(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, Communicator& comm, Domain_2& dm) {
  int n_ghost = 0;
  comm.comm_before_cal_force(p_arr, cl, n_ghost);

  auto f1 = [](TNode *pi, TNode *pj) {
    pi->interact(*pj);
  };
  auto f2 = [&dm](TNode *pi, TNode *pj) {
    pi->interact(*pj, dm);
  };
  cl.for_each_pair_slow(f1, f2);
  comm.clear_padded_particles(cl, p_arr, n_ghost);
}

template <typename TNode, typename TFunc>
void cal_force(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, Communicator& comm, TFunc for_all_pair_force) {
  int n_ghost = 0;
  comm.comm_before_cal_force(p_arr, cl, n_ghost);
  for_all_pair_force();
  comm.clear_padded_particles(cl, p_arr, n_ghost);
}

template <typename TNode, typename UniFunc>
void integrate(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, UniFunc f_move, Communicator& comm) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    f_move(*it);
  }
  cl.recreate(p_arr);
  comm.comm_after_integration(p_arr, cl);
}

//template <typename TNode, typename UniFunc, typename T1, typename T2>
//void integrate(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, UniFunc f_move,
//               Communicator &comm, std::vector<T1> &n_arr, std::vector<Vec_2<T2>> &v_arr) {
//  const auto end = p_arr.end();
//  for (auto it = p_arr.begin(); it != end; ++it) {
//    f_move(*it);
//  }
//  cl.recreate(p_arr, n_arr, v_arr);
//  comm.comm_after_integration(p_arr, cl, n_arr, v_arr);
//}

template <typename TNode, typename TInteract, typename TIntegrate, typename TOut>
void run(std::vector<TNode> &p_arr, TInteract interact_all,
         TIntegrate integrate_all, TOut out, int n_step) {
  for (int i = 1; i <= n_step; i++) {
    interact_all(p_arr);

    integrate_all(p_arr);

    out(i, p_arr);
  }
}

//void run_birth_death(int gl_par_num, const Vec_2<double>& gl_l,
//                     double eta, double alpha, unsigned long long seed,
//                     int n_step, int snap_interval, int snap_block_size,
//                     int box_len_birth);

void run_test(int gl_par_num, const Vec_2<double>&gl_l, double eta,
              unsigned long long seed, int n_step, const char* file_in = nullptr);