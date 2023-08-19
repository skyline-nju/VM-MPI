#pragma once
#include "config.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "exporter2D.h"
#include <iomanip>
#include "communicator2D.h"

void set_particle_num(int gl_par_num, int& my_par_num, int& my_par_num_max, MPI_Comm group_comm);

template <typename TNode, typename TRan, typename TDomain>
void ini_rand(std::vector<TNode>& p_arr, int gl_par_num, TRan& myran,
              CellListNode_2<TNode>& cl, const TDomain &dm) {
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
}

template <typename TNode, typename TFunc>
void cal_force(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, Communicator_2& comm, TFunc for_all_pair_force) {
  int n_ghost = 0;
  comm.comm_before_cal_force(p_arr, cl, n_ghost);
  for_all_pair_force();
  comm.clear_padded_particles(cl, p_arr, n_ghost);
}


// recreate cell lists when all particle have moved forward one step
template <typename TNode, typename UniFunc>
void integrate(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, UniFunc f_move, Communicator_2& comm, bool thick_shell=false) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    f_move(*it);
  }
  cl.recreate(p_arr, thick_shell);

  comm.comm_after_integration(p_arr, cl, thick_shell);

}

// update cell list once one particle has moved from one cell to another cell
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

template <typename TNode, typename TDomain>
double get_band_location(const std::vector<TNode>& p_arr, const TDomain& dm) {
  double my_m[2] = {0., 0.};
  double Lx = dm.gl_l().x;
  int np = p_arr.size();
  for (int i = 0; i < np; i++) {
    double theta = p_arr[i].pos.x / Lx * PI * 2;
    my_m[0] += std::cos(theta);
    my_m[1] += std::sin(theta);
  }
  double tot_m[2] = {0., 0.};
  MPI_Reduce(my_m, tot_m, 2, MPI_DOUBLE, MPI_SUM, 0, dm.comm());
  int my_rank;
  MPI_Comm_rank(dm.comm(), &my_rank);

  double xm;
  if (my_rank == 0) {
    double theta_m = std::atan2(tot_m[1], tot_m[0]);
    xm = theta_m / (PI * 2) * Lx;
    if (xm < 0) {
      xm += Lx;
    } else if (xm >= Lx) {
      xm -= Lx;
    }
  }
  MPI_Bcast(&xm, 1, MPI_DOUBLE, 0, dm.comm());
  return xm;
}

template <typename TDomain>
double get_left_edge(double xm, double width, const TDomain& dm) {
  double Lx = dm.gl_l().x;
  double x_left = xm + Lx / 2 - width / 2;
  if (x_left >= Lx)
    x_left -= Lx;
  return x_left;
}

template <typename TNode, typename TDomain>
void tag_particles(std::vector<int>& par_idx, double xmin, double width, const std::vector<TNode>& p_arr, const TDomain& dm) {
  int np = p_arr.size();
  double Lx = dm.gl_l().x;
  double xmax = xmin + width;
  if (xmax >= Lx)
    xmax -= Lx;
  if (xmin < xmax) {
    for (int i = 0; i < np; i++) {
      double x = p_arr[i].pos.x;
      if (x > xmin && x <= xmax)
        par_idx.push_back(i);
    }
  } else {
    for (int i = 0; i < np; i++) {
      double x = p_arr[i].pos.x;
      if (x > xmin || x <= xmax)
        par_idx.push_back(i);
    }
  } 
}

template <typename TNode, typename TDomain, typename TRan>
void equili_tagged_par(std::vector<int>& par_idx, double xmin, double width, TRan& myran, std::vector<TNode>& p_arr,
                     const TDomain& dm, CellListNode_2<TNode>& cl) {
  int my_np = par_idx.size();
  int tot_np;
  MPI_Reduce(&my_np, &tot_np, 1, MPI_INT, MPI_SUM, 0, dm.comm());
  int my_rank, proc_size;
  MPI_Comm_rank(dm.comm(), &my_rank);
  MPI_Comm_size(dm.comm(), &proc_size);

  int* np_arr = new int[proc_size] {};
  if (my_rank == 0) {
    int n_mean = tot_np / proc_size;
    for (int i = 0; i < proc_size; i++) {
      np_arr[i] = n_mean;
    }
    int j = int(myran.doub() * proc_size);
    np_arr[j] = tot_np - n_mean * (proc_size - 1);
  }

  int new_np;
  MPI_Scatter(np_arr, 1, MPI_INT, &new_np, 1, MPI_INT, 0, dm.comm());
  //std::cout << "excess particle number: " << my_np - new_np << std::endl;

  if (my_np > new_np) {
    std::vector<int> vac_pos;
    int n_removed = my_np - new_np;
    vac_pos.reserve(n_removed);
    shuffle(par_idx, myran, n_removed);
    do {
      int ip = par_idx.back();
      par_idx.pop_back();
      vac_pos.push_back(ip);
      cl.remove_node(p_arr[ip]);
    } while (vac_pos.size() < n_removed);
    std::sort(vac_pos.begin(), vac_pos.end(), std::greater<int>());
    cl.make_compact(p_arr, vac_pos);
    //std::cout << my_np - new_np << " particles removed" << std::endl;
  } else if (my_np < new_np) {
    int np_old = p_arr.size();
    do {
      double theta = myran.doub() * PI * 2;
      Vec_2<double> ori(std::cos(theta), std::sin(theta));
      Vec_2<double> pos;
      pos.x = xmin + width * myran.doub();
      if (pos.x >= dm.gl_l().x) {
        pos.x -= dm.gl_l().x;
      }
      pos.y = dm.origin().y + dm.l().y * myran.doub();
      p_arr.emplace_back(pos, ori);
      cl.add_node(p_arr.back());
    } while (p_arr.size() - np_old < new_np - my_np);
    //std::cout << p_arr.size() - np_old << " particles added" << std::endl;
  }

  delete[] np_arr;
}

template <typename TNode, typename TDomain, typename TRan>
void uniformize_gas(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl, const TDomain& dm, TRan& myran, 
                    double x_left, double tagged_region_width, double rho_0) {
  double ly = dm.l().y;
  std::vector<int> tagged_par_idx;
  //double x_band = get_band_location(p_arr, dm);
  //double x_left = get_left_edge(x_band, tagged_region_width, dm);
  int n_res = tagged_region_width * ly * rho_0 * 4;
  tagged_par_idx.reserve(n_res);
  tag_particles(tagged_par_idx, x_left, tagged_region_width, p_arr, dm);
  equili_tagged_par(tagged_par_idx, x_left, tagged_region_width, myran, p_arr, dm, cl);
  //std::vector<int> tagged_par_idx2;
  //tagged_par_idx2.reserve(n_res);
  //tag_particles(tagged_par_idx2, x_left, tagged_region_width, p_arr, dm);

  ////std::cout << tagged_par_idx2.size() << std::endl;
  //int np = p_arr.size();
  //int tot_np;
  //MPI_Reduce(&np, &tot_np, 1, MPI_INT, MPI_SUM, 0, dm.comm());
  //int my_rank;
  //MPI_Comm_rank(dm.comm(), &my_rank);
  //if (my_rank == 0)
  //  std::cout << "tot np: " << tot_np << ", band loc: " << x_band << std::endl;

}


template <typename TPar, typename TRan, typename TDomain>
void ini_particles(int gl_par_num, std::vector<TPar>& p_arr, 
                   const std::string& ini_mode, TRan& myran, int seed2, 
                   CellListNode_2<TPar>& cl, const TDomain& dm, int & t_beg) {
  if (ini_mode == "rand" || ini_mode == "ordered") {
#ifdef DISORDER_ON
    Ranq1 myran2(seed2);
    ini_rand(p_arr, gl_par_num, myran2, cl, dm);
#else
    ini_rand(p_arr, gl_par_num, myran, cl, dm);
#endif
    if (ini_mode == "ordered") {
      double angle = seed2 / 180. * PI;
      Vec_2<double> ori0 = Vec_2<double>(cos(angle), sin(angle));
      for (auto& p : p_arr) {
        p.ori = ori0;
#ifndef CONTINUE_DYNAMIC
        p.ori_next = ori0;
#else
        p.tau = 0.;
#endif
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

void run(int gl_par_num, const Vec_2<double>& gl_l,
         double D, double h, double v0,
         int n_step, std::string& ini_mode,
         unsigned long long seed,
         int snap_dt, int field_dt, int field_dx,
         MPI_Comm group_comm, MPI_Comm root_comm=MPI_COMM_WORLD);
