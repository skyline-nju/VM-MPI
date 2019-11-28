#pragma once
#include "config.h"
#include <iostream>
#include <vector>
#include "comn.h"
#include "domain2D.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

template <typename TRan>
void set_random_torque(double* theta, int n, double epsilon, TRan& myran) {
  const double d = 1.0 / (n - 1);
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

  template <typename Par>
  double get_torque(const Par& p);

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

  delete[] n_grids_v;
  delete[] displs;
  delete[] gl_theta;
#else
  set_random_torque(torque_, n_grids, epsilon, myran);
#endif
}

template <typename Par>
double RandTorque_2::get_torque(const Par& p) {
  const Vec_2<double> r = p.pos - origin_;
  return torque_[int(r.x) + int(r.y) * nx_];
}
