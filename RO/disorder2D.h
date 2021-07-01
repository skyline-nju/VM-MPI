#pragma once
#include "config.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "comn.h"
#include "domain2D.h"
#include "cellList2D.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

template <typename TRan>
void set_random_torque(double* theta, int n, double epsilon, TRan& myran) {
#ifdef RANDOM_FIELD
    const double d = 1.0 / n;
#else
    const double d = 1.0 / (n - 1);
#endif
  for (int i = 0; i < n; i++) {
    theta[i] = (-0.5 + i * d) * 2.0 * PI * epsilon;
  }
  shuffle(theta, n, myran);
}

class RandTorque_2 {
public:
  template <typename TRan>
  RandTorque_2(const double epsilon, TRan& myran, const Grid_2& grid, MPI_Comm group_comm);

  ~RandTorque_2() { delete[] torque_; }

  template <typename TPar>
  double get_torque(const TPar& p) const {
    const Vec_2<double> r = p.pos - origin_;
    return torque_[int(r.x) + int(r.y) * nx_];
  }

  double get_torque(int idx) const { return torque_[idx]; }
private:
  double* torque_;
  Vec_2<double> origin_;  // origin of the (sub)domain
  int nx_; // number of columns in x direction
};

template<typename TRan>
RandTorque_2::RandTorque_2(const double epsilon, TRan& myran, const Grid_2& grid, MPI_Comm group_comm)
                          : origin_(grid.lc() * grid.origin()), nx_(grid.n().x) {
  const int n_grids = nx_ * grid.n().y;
  torque_ = new double[n_grids];
#ifndef USE_MPI
  set_random_torque(torque_, n_grids, epsilon, myran);
#else
  const int gl_n_grids = grid.gl_n().x * grid.gl_n().y;
  double* gl_theta = nullptr;
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
  if (my_rank == 0) {
    gl_theta = new double[gl_n_grids];
    set_random_torque(gl_theta, gl_n_grids, epsilon, myran);
    double sum = 0;
    for (int i = 0; i < gl_n_grids; i++) {
      sum += gl_theta[i];
    }
    std::cout << "sum of torque = " << sum << std::endl;
  }

  //MPI_Scatter(gl_theta, n_grids, MPI_DOUBLE, torque_, n_grids, MPI_DOUBLE, 0, group_comm);
  //MPI_Barrier(group_comm);

  int *n_grids_v = new int[tot_proc];
  int *displs = new int[tot_proc];
  for (int i = 0; i < tot_proc; i++) {
    displs[i] = 0;
  }

  MPI_Gather(&n_grids, 1, MPI_INT, n_grids_v, 1, MPI_INT, 0, group_comm);
  if (my_rank == 0) {
    for (int i = 1; i < tot_proc; i++) {
      displs[i] = displs[i - 1] + n_grids_v[i - 1];
    }
  }
  MPI_Scatterv(gl_theta, n_grids_v, displs, MPI_DOUBLE, torque_, n_grids,
               MPI_DOUBLE, 0, group_comm);
  delete[] n_grids_v;
  delete[] displs;


  { // test
    double sum = 0;
    for (int i = 0; i < n_grids; i++) {
      sum += torque_[i];
    }
    double gl_sum = 0;
    MPI_Reduce(&sum, &gl_sum, 1, MPI_DOUBLE, MPI_SUM, 0, group_comm);
    if (my_rank == 0) {
      std::cout << "sum of torque = " << gl_sum << std::endl;
    }
  }

  delete[] gl_theta;
#endif
}

class RandField_2 {
public:
  template <typename TRan>

  RandField_2(double epsilon, TRan& myran, const Grid_2& grid, MPI_Comm group_comm);

  ~RandField_2() { delete[] field_; }

  template <typename TPar>
  const Vec_2<double>& get_field(const TPar& p) const {
    const Vec_2<double> r = p.pos - origin_;
    return field_[int(r.x) + int(r.y) * nx_];
  }

  template <typename TPar>
  void apply_field(TPar& p) const {
    const Vec_2<double> r = p.pos - origin_;
    const int idx = int(r.x) + int(r.y) * nx_;
    p.ori_next.x += p.n_neighb * field_[idx].x;
    p.ori_next.y += p.n_neighb * field_[idx].y;
  }
private:
  Vec_2<double>* field_;
  Vec_2<double> origin_;  // origin of the (sub)domain
  int nx_; // number of columns in x direction
};

template<typename TRan>
RandField_2::RandField_2(double epsilon, TRan& myran, const Grid_2& grid, MPI_Comm group_comm)
  : origin_(grid.lc()* grid.origin()), nx_(grid.n().x) {
  const int n_grids = nx_ * grid.n().y;
  field_ = new Vec_2<double>[n_grids];
  RandTorque_2 rand_torque(1., myran, grid, group_comm);
  double v_sum[2];
  v_sum[0] = v_sum[1] = 0.;
  for (int i = 0; i < n_grids; i++) {
    double theta = rand_torque.get_torque(i);
    field_[i].x = epsilon * cos(theta);
    field_[i].y = epsilon * sin(theta);
    v_sum[0] += field_[i].x;
    v_sum[1] += field_[i].y;
  }
#ifndef USE_MPI
  std::cout << "sum of random fields: " << v_sum[0] << "\t" << v_sum[1] << std::endl;
#else
  double v_sum_gl[2];
  MPI_Reduce(v_sum, v_sum_gl, 2, MPI_DOUBLE, MPI_SUM, 0, group_comm);
  int my_proc;
  MPI_Comm_rank(group_comm, &my_proc);
  if (my_proc == 0) {
    std::cout << "sum of random fields: " << v_sum_gl[0] << "\t" << v_sum_gl[1] << std::endl;
  }
#endif
}

class DiluteScatter: public CellListBase_2 {
public:
  template <class TDomain, class TGrid, class TRan>
  DiluteScatter(const TDomain& dm, const TGrid& grid,
                int n_ob, TRan& myran,
                MPI_Comm group_comm, const char* outfile);
  
  template <typename T>
  void add_obstacle(T x, T y);

  template <typename TPar>
  bool within(const TPar& p) const;

#ifndef DILUTE_COUPLING
  template <typename TPar>
#ifndef CONTINUE_DYNAMIC
  void scattering(TPar& p, double eps) const;
#else
  double scattering(const TPar& p) const;
#endif
#endif

private:
  Vec_2<double> half_gl_l_;
  std::vector<std::vector<Vec_2<double>>> obstacles_;
};

template <class TDomain, class TGrid, class TRan>
DiluteScatter::DiluteScatter(const TDomain& dm, const TGrid& grid,
                             int n_ob, TRan& myran,
                             MPI_Comm group_comn,
                             const char* outfile)
  : CellListBase_2(dm, grid, 1), half_gl_l_(gl_l_ * 0.5), obstacles_(n_.x * n_.y) {
  double *buf = new double[n_ob * 2];
  int my_rank = 0;
#ifdef USE_MPI
  MPI_Comm_rank(group_comn, &my_rank);
#endif
  if (my_rank == 0) {
    for (int i = 0; i < n_ob; i++) {
      buf[2 * i] = gl_l_.x * myran.doub();
      buf[2 * i + 1] = gl_l_.y * myran.doub();
    }
    std::ofstream fout(outfile, std::ios::binary);
    fout.write((char*)buf, sizeof(double) * n_ob * 2);
    fout.close();
  }
#ifdef USE_MPI
  MPI_Bcast(buf, n_ob * 2, MPI_DOUBLE, 0, group_comn);
  double xmin = origin_.x;
  double xmax = origin_.x + l_.x;
  double ymin = origin_.y;
  double ymax = origin_.y + l_.y;
#endif
  for (int i = 0; i < n_ob; i++) {
    double x = buf[2 * i];
    double y = buf[2 * i + 1];
#ifndef USE_MPI
    add_obstacle(x, y);
#else
    if (x - gl_l_.x >= xmin) {
      x -= gl_l_.x;
    } else if (x + gl_l_.x < xmax) {
      x += gl_l_.x;
    }
    if (y - gl_l_.y >= ymin) {
      y -= gl_l_.y;
    } else if (y + gl_l_.y < ymax) {
      y += gl_l_.y;
    }
    if (x >= xmin && x < xmax && y >= ymin && y < ymax) {
      add_obstacle(x, y);
    }
#endif
  }
  for (int ic = 0; ic < n_.x * n_.y; ic++) {
    obstacles_[ic].shrink_to_fit();
  }
  delete [] buf;
}


template <typename T>
void DiluteScatter::add_obstacle(T x, T y) {
  int col_c = get_nx(x);
  int row_c = get_ny(y);
  for (int drow = -1; drow < 2; drow++) {
    double dy = 0;
    int row_new = row_c + drow;
    if (!flag_ext_.y) {
      if (row_new < 0) {
        row_new += n_.y;
        dy = gl_l_.y;
      } else if (row_new >= n_.y) {
        row_new -= n_.y;
        dy = -gl_l_.y;
      }
    } else if (row_new < 0 || row_new >= n_.y) {
      continue;
    }
    int nx_row = n_.x * row_new;
    for (int dcol = -1; dcol < 2; dcol++) {
      double dx = 0;
      int col_new = col_c + dcol;
      if (!flag_ext_.x) {
        if (col_new < 0) {
          col_new += n_.x;
          dx = gl_l_.x;
        } else if (col_new >= n_.x) {
          col_new -= n_.x;
          dx = -gl_l_.x;
        }
      } else if (col_new < 0 || col_new >= n_.x) {
        continue;
      }
      int ic = col_new + nx_row;
      obstacles_[ic].push_back(Vec_2<double>(x+dx, y+dy));
    }
  }
}

template<typename TPar>
bool DiluteScatter::within(const TPar& p) const {
  int ic = get_ic(p);
  int max_ob = obstacles_[ic].size();
  bool ret = false;
  for (int io = 0; io < max_ob; io++) {
    double dx = p.pos.x - obstacles_[ic][io].x;
    double dy = p.pos.y - obstacles_[ic][io].y;
    double dis2 = dx * dx + dy * dy;
    if (dis2 < 1.) {
      ret = true;
      break;
    }
  }
  return ret;
}

#ifndef DILUTE_COUPLING

template <typename TPar>
#ifndef CONTINUE_DYNAMIC
void DiluteScatter::scattering(TPar &p, double eps) const{
  Vec_2<double> scat_vec(0., 0.);
  int ic = get_ic(p);
  int n_ob = 0;
  int max_ob = obstacles_[ic].size();
  for (int io=0; io < max_ob; io++) {
    double dx = p.pos.x - obstacles_[ic][io].x;
    double dy = p.pos.y - obstacles_[ic][io].y;
    double dis2 = dx * dx + dy * dy;
    if (dis2 < 1.) {
      double inv_dis = 1. / sqrt(dis2);
      scat_vec.x += dx * inv_dis;
      scat_vec.y += dy * inv_dis;
      n_ob++;
    }
  }
  if (n_ob > 0) {

#ifdef SCATTER_NORMED
      p.ori_next.x = n_ob * p.ori_next.x + eps * p.n_neighb * scat_vec.x;
      p.ori_next.y = n_ob * p.ori_next.y + eps * p.n_neighb * scat_vec.y;
#else
      p.ori_next.x = p.ori_next.x + eps * p.n_neighb * scat_vec.x;
      p.ori_next.y = p.ori_next.y + eps * p.n_neighb * scat_vec.y;
#endif
  }
  p.n_neighb = 1;
}
#else
double DiluteScatter::scattering(const TPar& p) const {
  double torque = 0.;
  int ic = get_ic(p);
  int n_ob = 0;
  int max_ob = obstacles_[ic].size();
  for (int io = 0; io < max_ob; io++) {
    double dx = p.pos.x - obstacles_[ic][io].x;
    double dy = p.pos.y - obstacles_[ic][io].y;
    double dis2 = dx * dx + dy * dy;
    if (dis2 < 1.) {
      double inv_dis = 1. / sqrt(dis2);
      dx *= inv_dis;
      dy *= inv_dis;
      torque += dy * p.ori.x - dx * p.ori.y;
      n_ob++;
    }
  }
  if (n_ob > 0) {
    torque /= n_ob;
  }
  return torque;
}
#endif
#endif