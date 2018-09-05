#pragma once
#define CHECK_RANGE
#include <vector>
#include "comn.h"
#include "mpi.h"

class RandTorqueArr_2 {
public:
  template <typename TRan>
  RandTorqueArr_2(double epsilon, TRan &myran, const Vec_2<double> &origin,
    const Vec_2<int> &cells_size, const Vec_2<int> &gl_cells_size);

  ~RandTorqueArr_2() { delete[] torque_; }

  template <typename TPar>
  double get_torque(const TPar &par) const;
private:
  double *torque_;
  Vec_2<double> origin_;
  int nx_;
  int n_tot_;
};

template<typename TRan>
RandTorqueArr_2::RandTorqueArr_2(double epsilon, TRan & myran,
                                 const Vec_2<double>& origin,
                                 const Vec_2<int>& cells_size,
                                 const Vec_2<int>& gl_cells_size)
  : origin_(origin) {
  nx_ = cells_size.x;
  n_tot_ = nx_ * cells_size.y;
  int gl_n_tot = gl_cells_size.x * gl_cells_size.y;

  double *gl_torque = nullptr;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0) {
    double d = 1.0 / (gl_n_tot - 1);
    gl_torque = new double[gl_n_tot];
    for (int i = 0; i < gl_n_tot; i++) {
      gl_torque[i] = (-0.5 + i * d) * epsilon * 2.0 * PI;
    }
    shuffle(gl_torque, gl_n_tot, myran);
  }

  torque_ = new double[n_tot_];
  MPI_Scatter(gl_torque, n_tot_, MPI_DOUBLE,
              torque_, n_tot_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  delete[] gl_torque;

  double torque_sum = 0;
  for (int i = 0; i < n_tot_; i++) {
    torque_sum += torque_[i];
  }
  double gl_torque_sum;
  MPI_Reduce(&torque_sum, &gl_torque_sum, 1, MPI_DOUBLE,
             MPI_SUM, 0, MPI_COMM_WORLD);
  if (my_rank == 0) {
    std::cout << "sum of random torque = " << gl_torque_sum << std::endl;
  }
}

template <typename TPar>
double RandTorqueArr_2::get_torque(const TPar& par) const {
  const Vec_2<double> r = par.pos - origin_;
  const int idx = int(r.x) + int(r.y) * nx_;
#ifdef CHECK_RANGE
  if (idx < 0 || idx > n_tot_) {
    std::cerr << "out of range when getting torque\n";
    std::cerr << "idx = " << idx << std::endl;
    std::cerr << "origin = " << origin_ << std::endl;
    std::cerr << "pos = " << par.pos << std::endl;
    exit(-1);
  }
#endif
  return torque_[idx];
}
