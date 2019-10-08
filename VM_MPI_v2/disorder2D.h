#pragma once
#include "config.h"
#include <vector>
#include "comn.h"
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
               const Vec_2<int>& cells_size, const Vec_2<int>& gl_cells_size);

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
                           const Vec_2<int>& cells_size, const Vec_2<int>& gl_cells_size)
                           :origin_(origin), nx_(cells_size.x) {
  const int n_tot = nx_ * cells_size.y;
  torque_ = new double[n_tot];
#ifdef USE_MPI
  const int gl_n_tot = gl_cells_size.x * gl_cells_size.y;
  double* gl_theta = nullptr;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0) {
    gl_theta = new double[gl_n_tot];
    set_random_torque(gl_theta, gl_n_tot, epsilon, myran);
  }
  MPI_Scatter(gl_theta, n_tot, MPI_DOUBLE, torque_,
              n_tot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  delete[] gl_theta;
#else
  set_random_torque(torque_, n_tot, epsilon, myran);
#endif
}

template <typename Par>
double RandTorque_2::get_torque(const Par& p) {
  const Vec_2<double> r = p.pos - origin_;
  return torque_[int(r.x) + int(r.y) * nx_];
}
