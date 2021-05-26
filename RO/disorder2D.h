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
  RandTorque_2(const double epsilon, TRan& myran, const Vec_2<double>& origin,
               const Vec_2<int>& cells_size, const Vec_2<int>& gl_cells_size, MPI_Comm group_comm);

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
RandTorque_2::RandTorque_2(const double epsilon, TRan& myran, const Vec_2<double>& origin,
                           const Vec_2<int>& cells_size, const Vec_2<int>& gl_cells_size,
                           MPI_Comm group_comm)
                           :origin_(origin), nx_(cells_size.x) {
  const int n_tot = nx_ * cells_size.y;
  torque_ = new double[n_tot];
#ifdef USE_MPI
  const int gl_n_tot = gl_cells_size.x * gl_cells_size.y;
  double* gl_theta = nullptr;
  int my_rank;
  MPI_Comm_rank(group_comm, &my_rank);
  if (my_rank == 0) {
    gl_theta = new double[gl_n_tot];
    set_random_torque(gl_theta, gl_n_tot, epsilon, myran);
    double sum = 0;
    for (int i = 0; i < gl_n_tot; i++) {
      sum += gl_theta[i];
    }
    std::cout << "sum of torque = " << sum << std::endl;
  }

  MPI_Scatter(gl_theta, n_tot, MPI_DOUBLE, torque_,
              n_tot, MPI_DOUBLE, 0, group_comm);
  { // test
    double sum = 0;
    for (int i = 0; i < n_tot; i++) {
      sum += torque_[i];
    }
    double gl_sum = 0;
    MPI_Reduce(&sum, &gl_sum, 1, MPI_DOUBLE, MPI_SUM, 0, group_comm);
    if (my_rank == 0) {
      std::cout << "sum of torque = " << gl_sum << std::endl;
    }
  }
  delete[] gl_theta;
#else
  set_random_torque(torque_, n_tot, epsilon, myran);
#endif
}

template<typename TRan>
RandTorque_2::RandTorque_2(const double epsilon, TRan& myran, const Grid_2& grid, MPI_Comm group_comm)
                          : origin_(grid.lc() * grid.origin()), nx_(grid.n().x) {
  const int n_grids = nx_ * grid.n().y;
  torque_ = new double[n_grids];
#ifdef USE_MPI
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
#else
  set_random_torque(torque_, n_grids, epsilon, myran);
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

  double v_sum_gl[2];
  MPI_Reduce(v_sum, v_sum_gl, 2, MPI_DOUBLE, MPI_SUM, 0, group_comm);
  int my_proc;
  MPI_Comm_rank(group_comm, &my_proc);
  if (my_proc == 0) {
    std::cout << "sum of random fields: " << v_sum_gl[0] << "\t" << v_sum_gl[1] << std::endl;
  }
}

class RandScatter: public CellListBase_2 {
public:
  template <class TDomain, class TGrid, class TRan>
  RandScatter(const TDomain& dm, const TGrid& grid, int n_ob, TRan& myran, MPI_Comm group_comm, const char* outfile);

  template <typename T>
  void scattering(T x, T y, Vec_2<double>& scat_vec, int &n_ob) const;

  template <typename TPar>
  void scattering(TPar& p, double eps) const;
  
  int get_col(int col0, int dcol) const;
  int get_row(int row0, int drow) const;

private:
  Vec_2<double> half_gl_l_;
  std::vector<std::vector<Vec_2<double>>> obstacles_;
};

template <class TDomain, class TGrid, class TRan>
RandScatter::RandScatter(const TDomain& dm, const TGrid& grid, int n_ob, TRan& myran, MPI_Comm group_comm, const char* outfile)
    : CellListBase_2(dm, grid, 1), half_gl_l_(gl_l_ * 0.5), obstacles_(n_.x * n_.y) {
  int my_rank;
  MPI_Comm_rank(group_comm, &my_rank);
  double *buf = new double[n_ob * 2];
  if (my_rank == 0) {
    for (int i = 0; i < n_ob; i++) {
      buf[2 * i] = gl_l_.x * myran.doub();
      buf[2 * i + 1] = gl_l_.y * myran.doub();
    }
    std::ofstream fout(outfile, std::ios::binary);
    fout.write((char*)buf, sizeof(double) * n_ob * 2);
    fout.close();
  }
  MPI_Bcast(buf, n_ob * 2, MPI_DOUBLE, 0, group_comm);

  double xmin = origin_.x;
  double xmax = origin_.x + l_.x;
  double ymin = origin_.y;
  double ymax = origin_.y + l_.y;
  for (int i = 0; i < n_ob; i++) {
    double x = buf[2 * i];
    double y = buf[2 * i + 1];
    if (x >= xmin && x < xmax && y >= ymin && y < ymax) {
      int ic = get_nx(x) + n_.x * get_ny(y);
      obstacles_[ic].push_back(Vec_2<double>(x, y));
    }
  }
  delete [] buf;
}

inline int RandScatter::get_col(int col0, int dcol) const {
  int col = col0 + dcol;
  if (!flag_ext_.x) {
    if (col < 0) {
      col += n_.x;
    } else if (col >= n_.x) {
      col -= n_.x;
    }
  }
  return col;
}
inline int RandScatter::get_row(int row0, int drow) const {
  int row = row0 + drow;
  if (!flag_ext_.y) {
    if (row < 0) {
      row += n_.y;
    } else if (row >= n_.y) {
      row -= n_.y;
    }
  }
  return row;
}

template <typename T>
void RandScatter::scattering(T x, T y, Vec_2<double>& scat_vec, int &n_ob) const {
  int col0 = get_nx(x);
  int row0 = get_ny(y);
  for (int drow = -1; drow < 2; drow++) {
    int row = get_row(row0, drow);
    int row_nx = row * n_.x;
    for (int dcol = -1; dcol < 2; dcol++) {
      int col = get_col(col0, dcol);
      int ic = col + row_nx;
      int size_ob = obstacles_[ic].size();
      for (int io = 0; io < size_ob; io++) {
        double dx = x - obstacles_[ic][io].x;
        double dy = y - obstacles_[ic][io].y;
        if (!flag_ext_.x) {
          untangle_1(dx, gl_l_.x, half_gl_l_.x);
        }
        if (!flag_ext_.y) {
          untangle_1(dy, gl_l_.y, half_gl_l_.y);
        }
        double dis2 = dx * dx + dy * dy;
        if (dis2 < 1.) {
          double inv_dis = 1. / sqrt(dis2);
          scat_vec.x += dx * inv_dis;
          scat_vec.y += dy * inv_dis;
          n_ob++;
        }
      }
    }
  }
}

template <typename TPar>
void RandScatter::scattering(TPar& p, double eps) const {
  int n_ob = 0;
  Vec_2<double> scat_vec(0., 0.);
  scattering(p.pos.x, p.pos.y, scat_vec, n_ob);
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