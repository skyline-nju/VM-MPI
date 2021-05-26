#pragma once
#include "config.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "exporter2D.h"
#include <iomanip>
#include "communicator2D.h"

int ini_my_par_num(int gl_par_num, MPI_Comm group_comm);

void set_particle_num(int gl_par_num, int& my_par_num, int& my_par_num_max, MPI_Comm group_comm);

template <typename TNode, typename TRan, typename TDomain>
void ini_rand(std::vector<TNode>& p_arr, int gl_par_num, TRan& myran,
              CellListNode_2<TNode>& cl, TDomain &dm) {
  int my_par_num, n_max;
  set_particle_num(gl_par_num, my_par_num, n_max, dm.comm());
  p_arr.reserve(n_max);
  for (int i = 0; i < my_par_num; i++) {
    p_arr.emplace_back(myran, dm.l(), dm.origin());
  }
  cl.create(p_arr);
}

template <typename TNode, typename TDomain>
void ini_from_snap(std::vector<TNode>& p_arr, int gl_par_num,
                   CellListNode_2<TNode>& cl,
                   TDomain &dm, const char* filename) {
  int my_rank;
  MPI_Comm_rank(dm.comm(), &my_rank);

  int my_par_num, n_max;
  set_particle_num(gl_par_num, my_par_num, n_max, dm.comm());

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

template <typename TNode, typename TFunc>
void cal_force(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, Communicator_2& comm, TFunc for_all_pair_force) {
  int n_ghost = 0;
  comm.comm_before_cal_force(p_arr, cl, n_ghost);
  for_all_pair_force();
  comm.clear_padded_particles(cl, p_arr, n_ghost);
}

template <typename TNode, typename UniFunc>
void integrate(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, UniFunc f_move, Communicator_2& comm) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    f_move(*it);
  }
  cl.recreate(p_arr);

  comm.comm_after_integration(p_arr, cl);

}

template <typename TNode, typename UniFunc>
void integrate2(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, UniFunc f_move, Communicator_2& comm) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic_old = cl.get_ic(*it);
    f_move(*it);
    int ic_new = cl.get_ic(*it);
    if (ic_old != ic_new) {
      cl.update(*it, ic_old, ic_new);
    }
  }
  comm.comm_after_integration(p_arr, cl);
}

#if defined RANDOM_TORQUE || defined RANDOM_FIELD
void run_quenched(int gl_par_num, const Vec_2<double>& gl_l, double eta, double eps,
                  unsigned long long seed, int seed2, int n_step,
                  MPI_Comm group_comm, MPI_Comm root_comm = MPI_COMM_WORLD);
#endif


void run_RO(int gl_par_num, const Vec_2<double>& gl_l,
            double eta, double rho_s, double eps,
            unsigned long long seed, int seed2, int n_step,
            MPI_Comm group_comm, MPI_Comm root_comm);