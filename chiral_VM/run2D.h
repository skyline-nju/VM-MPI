#pragma once
#include "config.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "exporter2D.h"
#include <iomanip>
#ifdef USE_MPI
#include "communicator2D.h"
#endif

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
                   CellListNode_2<TNode>& cl, const TDomain& dm,
                   const char* filename) {
  int my_par_num, n_max;
  set_particle_num(gl_par_num, my_par_num, n_max, dm.comm());
  p_arr.reserve(n_max);
  float buf[3];
  std::ifstream fin(filename, std::ios::binary);
  double Lx = dm.gl_l().x;
  double Ly = dm.gl_l().y;
  for (int i = 0; i < gl_par_num; i++) {
    fin.read((char*)buf, sizeof(float) * 3);
    TNode p(buf);
    if (p.pos.x < 0) {
      p.pos.x += Lx;
    } else if (p.pos.x >= Lx) {
      p.pos.x -= Lx;
    }
    if (p.pos.y < 0) {
      p.pos.y += Ly;
    } else if (p.pos.y >= Ly) {
      p.pos.y -= Ly;
    }
    if (dm.contain_particle(p)) {
      p_arr.push_back(p);
    }
  }
  fin.close();
  cl.create(p_arr);

#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(dm.comm(), &my_rank);
  int my_par = p_arr.size();
  int tot_par;
  MPI_Reduce(&my_par, &tot_par, 1, MPI_INT, MPI_SUM, 0, dm.comm());

  if (my_rank == 0) {
    if (tot_par != gl_par_num) {
      std::cout << "Error when loading " << filename << ", " << tot_par
        << " instead of " << gl_par_num << " have been loaded" << std::endl;
      exit(1);
    } else {
      std::cout << "Load " << tot_par << " from " << filename << std::endl;
    }
  }
#endif
}

#ifdef USE_MPI
template <typename TNode, typename TFunc>
void cal_force(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, Communicator_2& comm, TFunc for_all_pair_force) {
  int n_ghost = 0;
  comm.comm_before_cal_force(p_arr, cl, n_ghost);
  for_all_pair_force();
  comm.clear_padded_particles(cl, p_arr, n_ghost);
}
#endif


// recreate cell lists when all particle have moved forward one step
template <typename TNode, typename UniFunc>
#ifdef USE_MPI
void integrate(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, UniFunc f_move, Communicator_2& comm) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    f_move(*it);
  }
  cl.recreate(p_arr);

  comm.comm_after_integration(p_arr, cl);
}
#else
void integrate(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, UniFunc f_move) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    f_move(*it);
  }
  cl.recreate(p_arr);
}
#endif

// update cell list once one particle has moved from one cell to another cell
template <typename TNode, typename UniFunc>
#ifdef USE_MPI
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
#else
void integrate2(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, UniFunc f_move) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic_old = cl.get_ic(*it);
    f_move(*it);
    int ic_new = cl.get_ic(*it);
    if (ic_old != ic_new) {
      cl.update(*it, ic_old, ic_new);
    }
  }
}
#endif

template <typename TPar, typename TRan, typename TDomain>
void ini_particles(int gl_par_num, std::vector<TPar>& p_arr,
                   const std::string& ini_mode, TRan& myran,
                   CellListNode_2<TPar>& cl, const TDomain& dm, int& t_beg) {
  if (ini_mode == "rand" || ini_mode == "ordered") {
    unsigned long long my_seed = myran.int64();
    Ranq1 myran2(my_seed);
    ini_rand(p_arr, gl_par_num, myran2, cl, dm);
    if (ini_mode == "ordered") {
      Vec_2<double> ori0 = Vec_2<double>(1, 0);
      for (auto& p : p_arr) {
        p.ori = ori0;
        p.tau = 0.;
      }
    }
    t_beg = 0;
  } else {
    ini_from_snap(p_arr, gl_par_num, cl, dm, ini_mode.c_str());
    std::vector<std::string> str_vec = split(ini_mode, "_");
    std::vector<std::string> str_vec2 = split(str_vec.back(), ".");
    str_to_num(str_vec2[0], t_beg);
  }
}

template <typename UniFunc>
void run(UniFunc one_step, int n_step,
         MPI_Comm group_comm,
         MPI_Comm root_comm) {
#ifndef USE_MPI
  for (int i = 1; i <= n_step; i++) {
    one_step(i);
  }
#else
if (root_comm == MPI_COMM_WORLD) {
  for (int i = 1; i <= n_step; i++) {
    one_step(i);
  }
} else {
  int i_step = 1;
  MPI_Win win_root;
  int finished_group = 0;
  int n_group;
  int my_rank;
  MPI_Comm_rank(group_comm, &my_rank);
  if (my_rank == 0) {
    MPI_Win_create(&i_step, 1, sizeof(int), MPI_INFO_NULL, root_comm, &win_root);
    MPI_Comm_size(root_comm, &n_group);
    std::cout << "number of groups: " << n_group << std::endl;
  }
  MPI_Bcast(&n_group, 1, MPI_INT, 0, group_comm);
  while (finished_group < n_group) {
    one_step(i_step);
    i_step++;
    if (i_step > n_step && i_step % 100 == 1) {
      if (my_rank == 0) {
        finished_group = 0;
        for (int j = 0; j < n_group; j++) {
          int remote_i_step;
          MPI_Win_lock(MPI_LOCK_SHARED, j, 0, win_root);
          MPI_Get(&remote_i_step, 1, MPI_INT, j, 0, 1, MPI_INT, win_root);
          MPI_Win_unlock(j, win_root);
          if (remote_i_step > n_step) {
            finished_group++;
          }
        }
      }
      MPI_Bcast(&finished_group, 1, MPI_INT, 0, group_comm);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif
}