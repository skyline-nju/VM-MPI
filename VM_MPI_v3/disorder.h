/**
 * @brief  Quenched disorder
 * 
 * @file disorder.h
 * @author skyline-nju
 * @date 2018-08-30
 */

#pragma once

#include <vector>
#include "vect.h"
#include "comn.h"
#include "mpi.h"

class RandTorque {
public:
  RandTorque(): cos_theta(0), sin_theta(0), ori(Vec_3<double>()) {}

  RandTorque(double c, double s, const Vec_3<double> &v)
    : cos_theta(c), sin_theta(s), ori(v) {}

  template <typename TRan>
  RandTorque(double theta, TRan &myran);

  template <typename T>
  void rotate(Vec_3<T> &v) const;

  template <typename TRan>
  static void ini(double epsilon, TRan &myran, std::vector<RandTorque> &torque,
                  int size, int gl_size);

  double cos_theta;
  double sin_theta;
  Vec_3<double> ori;
};



template <typename TRan>
RandTorque::RandTorque(double theta, TRan &myran) {
  cos_theta = std::cos(theta);
  sin_theta = std::sin(theta);
  sphere_point_picking(ori.x, ori.y, ori.z, myran);
}

template <typename T>
void RandTorque::rotate(Vec_3<T>& v) const {
  v.rotate(cos_theta, sin_theta, ori);
}

template <typename TRan>
void RandTorque::ini(double epsilon, TRan& myran, std::vector<RandTorque>& torque_arr,
                     int size, int gl_size) {
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  std::vector<Vec_3<double>> ori_arr;
  std::vector<double> theta_arr;
  ori_arr.reserve(size);
  theta_arr.reserve(size);
  torque_arr.reserve(size);
  Vec_3<double> v_sum{};

  for (int i = 0; i < size; i++) {
    theta_arr.push_back(epsilon * PI * 2 * (myran.doub() - 0.5));
    Vec_3<double> v{};
    sphere_point_picking(v.x, v.y, v.z, myran);
    ori_arr.push_back(v);
    v_sum += v * theta_arr.back();
  }

  Vec_3<double> gl_v{};
  MPI_Reduce(&v_sum.x, &gl_v, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gl_v.x, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  gl_v.x /= gl_size;
  gl_v.y /= gl_size;
  gl_v.z /= gl_size;

  v_sum = Vec_3<double>();
  for (int i = 0; i < size; i++) {
    Vec_3<double> v = theta_arr[i] * ori_arr[i] - gl_v;
    theta_arr[i] = v.module();
    ori_arr[i] = v / theta_arr[i];
    double c = std::cos(theta_arr[i]);
    double s = std::sin(theta_arr[i]);
    torque_arr.emplace_back(c, s, ori_arr[i]);
    v_sum += theta_arr[i] * ori_arr[i];
  }
  MPI_Reduce(&v_sum.x, &gl_v, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (my_rank == 0) {
    std::cout << "sum of random torques:" << gl_v << std::endl;
  }
}

class RandTorqueArr {
public:
  template <typename TRan>
  RandTorqueArr(double epsilon, TRan &myran, const Vec_3<double> &orgin,
                const Vec_3<int> &cells_size, const Vec_3<int> &gl_cells_size);

  template <typename TPar>
  const RandTorque& get_torque(const TPar &par) const;

private:
  std::vector<RandTorque> arr_;
  Vec_3<double> origin_;
  int nx_;
  int nx_ny_;
  int n_tot_;
};

template <typename TRan>
RandTorqueArr::RandTorqueArr(double epsilon, TRan& myran, const Vec_3<double> &origin,
                             const Vec_3<int>& cells_size, const Vec_3<int>& gl_cells_size)
  : origin_(origin) {
  nx_ = cells_size.x;
  nx_ny_ = nx_ * cells_size.y;
  n_tot_ = nx_ny_ * cells_size.z;
  int gl_n_tot = gl_cells_size.x * gl_cells_size.y * gl_cells_size.z;
  RandTorque::ini(epsilon, myran, arr_, n_tot_, gl_n_tot);
}

template <typename TPar>
const RandTorque& RandTorqueArr::get_torque(const TPar &par) const {
  const Vec_3<double> r = par.pos - origin_;
  const int idx = int(r.x) + int(r.y) * nx_ + int(r.z) * nx_ny_;
  if (idx < 0 || idx > n_tot_) {
    std::cerr << "idx = " << idx << std::endl;
    std::cerr << "origin = " << origin_ << std::endl;
    std::cerr << "pos = " << par.pos << std::endl;
    exit(-1);
  }
  return arr_[idx];
}