#pragma once
#include "domain3D.h"
#include "cellList3D.h"
#include <iomanip>
#include "disorder.h"

int ini_my_par_num(int gl_par_num);

template <class TPar>
void cal_order_para(const std::vector<TPar> &p_arr,double &phi, Vec_3<double> &v_mean) {
  Vec_3 <double> sum_v{};
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    sum_v += (*it).ori;
  }
  int tot_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  if (tot_proc > 1) {
    Vec_3<double> gl_sum_v{};
    int my_np = p_arr.size();
    int gl_np = 0;
    MPI_Reduce(&sum_v.x, &gl_sum_v.x, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&my_np, &gl_np, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0) {
      v_mean = gl_sum_v / gl_np;
      phi = v_mean.module();
      std::cout << std::setprecision(10) << "phi = " << phi << "\tori =\t" << v_mean << "\tN = " << gl_np << std::endl;
    }
  } else {
    v_mean = sum_v / p_arr.size();
    phi = v_mean.module();
    std::cout << std::setprecision(16) << std::setw(18) << "phi = " << phi << "\t" << "ori = " << v_mean << std::endl;
  }
}

template <typename TNode, typename TRan>
void ini_rand(std::vector<TNode>& p_arr, int gl_par_num, TRan& myran,
              CellListNode_3<TNode>& cl, Domain_3 &dm) {
  const int my_par_num = ini_my_par_num(gl_par_num);
#ifdef USE_MPI
  p_arr.reserve(my_par_num * 3);
  dm.set_max_buf_size(gl_par_num, 20.);
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
  auto f1 = [](TNode *pi, TNode *pj) {
    pi->interact(*pj);
  };
  auto f2 = [&dm](TNode *pi, TNode *pj) {
    pi->interact(*pj, dm);
  };
  auto f3 = [](TNode *pi, TNode *pj, const Vec_3<double> &offset) {
    pi->interact(*pj, offset);
  };

#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int n_ghost = 0;

  comm_par_before_interact(dm.neighbor, dm.inner_shell, dm.max_buf_size(),
    p_arr, cl, n_ghost);
#endif
  //Vec_3<int> beg1 = Vec_3<int>();
  //Vec_3<int> end1 = Vec_3<int>(cl.cells_size().x, cl.cells_size().y, 1);
  //Vec_3<int> beg2 = Vec_3<int>(0, 0, 1);
  //Vec_3<int> end2 = Vec_3<int>(cl.cells_size().x, cl.cells_size().y, cl.cells_size().z - 2);
  //Vec_3<int> beg3 = Vec_3<int>(0, 0, cl.cells_size().z - 2);
  //Vec_3<int> end3 = Vec_3<int>(cl.cells_size().x, cl.cells_size().y, cl.cells_size().z - 1);

  //cl.for_each_pair(f1, f2, beg2, end2);
  //cl.for_each_pair(f1, f2, beg1, end1);
  //cl.for_each_pair(f1, f2, beg3, end3);

  //cl.for_each_pair_slow(f2, beg2, end2);
  //cl.for_each_pair_slow(f2, beg1, end1);
  //cl.for_each_pair_slow(f2, beg3, end3);

  //cl.for_each_pair_fast(f1, f3, beg2, end2);
  //cl.for_each_pair_fast(f1, f3, beg1, end1);
  //cl.for_each_pair_fast(f1, f3, beg3, end3);
  cl.for_each_pair_fast(f1, f3);
  //cl.for_each_pair(f1, f2);
  //cl.for_each_pair_slow(f2);


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

template <typename TNode, typename TRan>
void integrate(std::vector<TNode>& p_arr, TRan &myran,
               CellListNode_3<TNode>& cl, const Domain_3 &dm,
               double eta, double v0, const RandTorqueArr &torques) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    (*it).move(eta, v0, torques.get_torque(*it), myran, dm);
  }
  cl.recreate(p_arr);
#ifdef USE_MPI
  comm_par_after_move(dm.neighbor, dm.outer_shell, dm.max_buf_size(), p_arr, cl);
#endif
}

template <typename TNode, typename TInteract, typename TIntegrate>
void run(std::vector<TNode> &p_arr, int gl_par_num, TInteract interact,
         TIntegrate integrate, int n_step, int sep) {
  int my_rank, tot_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  auto t1 = std::chrono::system_clock::now();
  double phi_sum = 0;
  int count = 0;
  for (int i = 1; i <= n_step; i++) {
    
    interact(p_arr);

    integrate(p_arr);

    if (i % sep == 0) {
      double phi = 0;
      Vec_3<double> v_mean{};
      cal_order_para(p_arr, phi, v_mean);
      if (i > 100000) {
        phi_sum += phi;
        count++;
      }
    }
  }
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  if (my_rank == 0) {
    if (count > 0)
      std::cout << "phi = " << phi_sum / count << std::endl;  
    std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
    std::cout << std::scientific <<"speed: " << n_step * gl_par_num / elapsed_time.count() / tot_proc
      << " particle step per second per core" << std::endl;
  }
}

template <typename TNode, typename TInteract, typename TIntegrate, typename TOut>
void run(std::vector<TNode> &p_arr, int gl_par_num, TInteract interact,
         TIntegrate integrate, TOut out, int n_step) {
  int my_rank, tot_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  auto t1 = std::chrono::system_clock::now();
  for (int i = 1; i <= n_step; i++) {
    if (my_rank == 0 && i % 1000 == 0) {
      auto t2 = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_time = t2 - t1;
      std::cout << "t = " << i << "\telapsed seconds = " << elapsed_time.count() << std::endl;
    }
    interact(p_arr);

    integrate(p_arr);

    out(i, p_arr);
  }
}