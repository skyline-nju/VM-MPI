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
#ifdef RANDOM_TORQUE
    const double d = 1.0 / (n - 1);
#elif defined RANDOM_FIELD
    const double d = 1.0 / n;
#elif defined RANDOM_POTENTIAL
    const double d = epsilon / (n - 1);
#endif
  for (int i = 0; i < n; i++) {
#ifndef RANDOM_POTENTIAL
    theta[i] = (-0.5 + i * d) * 2.0 * PI * epsilon;
#else
    theta[i] = i * d * 2.0 * PI;
#endif
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

class RandPotential_2 {
public:
  template <typename TRan>
  RandPotential_2(double eta, double eps, TRan& myran, const Grid_2& grid, MPI_Comm group_comm);

  ~RandPotential_2() { delete[] potential_; }

  template <typename TPar>
  double get_potential(const TPar& p) const {
    const Vec_2<double> r = p.pos - origin_;
    return potential_[int(r.x) + int(r.y) * nx_];
  }


private:
  double* potential_;
  Vec_2<double> origin_;   // origin of the (sub)domain
  int nx_;                // number of columns in x direction
};

template<typename TRan>
RandPotential_2::RandPotential_2(double eta, double eps, TRan& myran, 
                                 const Grid_2& grid, MPI_Comm group_comm)
  : origin_(grid.lc()* grid.origin()), nx_(grid.n().x) {
  const int n_grids = nx_ * grid.n().y;
  potential_ = new double[n_grids];
  RandTorque_2 rand_torque(eps, myran, grid, group_comm);
  double p_sum = 0.;
  for (int i = 0; i < n_grids; i++) {
    double theta = rand_torque.get_torque(i);
    potential_[i] = theta + eta * 2.0 * PI;
    p_sum += theta;
  }
  double p_sum_gl;
  MPI_Reduce(&p_sum, &p_sum_gl, 1, MPI_DOUBLE, MPI_SUM, 0, group_comm);
  int my_proc;
  MPI_Comm_rank(group_comm, &my_proc);
  if (my_proc == 0) {
    std::cout << "mean random potential: " << p_sum_gl/n_grids << std::endl;
  }
}

