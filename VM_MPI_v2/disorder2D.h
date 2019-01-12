#pragma once
#include "config.h"
#include <vector>
#include "comn.h"
#ifdef USE_MPI
#include "mpi.h"
#endif


class Disorder_2 {
public:
  template <typename TRan>
  Disorder_2(const double epsilon, TRan &myran, const Vec_2<double> &origin,
           const Vec_2<int> &cells_size, const Vec_2<int> &gl_cells_size);

  const Vec_2<double>& get_disorder(const Vec_2<double> &pos) const;

  template <typename TRan>
  void ini_random_angle(double *theta, int n, TRan &myran) const;

  void eval(const Vec_2<double> &pos, Vec_2<double> &v) const;

private:
  std::vector<Vec_2<double>> u_arr_;
  Vec_2<double> origin_;
  int nx_;
};

template <typename TRan>
Disorder_2::Disorder_2(const double epsilon, TRan &myran, const Vec_2<double> &origin,
                       const Vec_2<int> &cells_size, const Vec_2<int> &gl_cells_size):
                       origin_(origin), nx_(cells_size.x) {
  const int n_tot = nx_ * cells_size.y;
  double *theta = new double[n_tot];
#ifdef USE_MPI
  const int gl_n_tot = gl_cells_size.x * gl_cells_size.y;
  double *gl_theta = nullptr;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0) {
    gl_theta = new double[gl_n_tot];
    ini_random_angle(gl_theta, gl_n_tot, myran);
  }
  MPI_Scatter(gl_theta, n_tot, MPI_DOUBLE, theta,
              n_tot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  delete[] gl_theta;
#else
  ini_random_angle(theta, n_tot, myran);
#endif
  u_arr_.reserve(n_tot);
  for (int i = 0; i < n_tot; i++) {
#ifdef RANDOM_TORQUE
    const double phi = epsilon * theta[i];
    const Vec_2<double> u(std::cos(phi), std::sin(phi));
#else
    const Vec_2<double> u(epsilon * std::cos(theta[i]), epsilon * std::sin(theta[i]));
#endif
    u_arr_.push_back(u);
  }
  delete[] theta;
   
  double sum_v[2] = { 0., 0. };
  for (auto &i : u_arr_) {
    sum_v[0] += i.x;
    sum_v[1] += i.y;
  }
#ifdef USE_MPI
  double gl_sum_v[2] = { 0., 0. };
  MPI_Reduce(sum_v, gl_sum_v, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  sum_v[0] = gl_sum_v[0];
  sum_v[1] = gl_sum_v[1];
#endif
  if (my_rank == 0) {
    std::cout << "sum of the random";
#ifdef RANDOM_TORQUE
    std::cout << "torque: ";
#elif RANDOM_FIELD
    std::cout << "field: ";
#else
    std::cout << "stress: ";
#endif
    std::cout << sum_v[0] << "\t" << sum_v[1] << std::endl;
    
  }
}

template <typename TRan>
void Disorder_2::ini_random_angle(double* theta, int n, TRan &myran) const{
  const double d = 1.0 / (n - 1);
  for (int i = 0; i < n; i++) {
    theta[i] = (-0.5 + i * d) * 2.0 * PI;
  }
  shuffle(theta, n, myran);
}


inline const Vec_2<double>& Disorder_2::get_disorder(const Vec_2<double> &pos) const {
  const Vec_2<double> r = pos - origin_;
  return u_arr_[int(r.x) + int(r.y) * nx_];
}

inline void Disorder_2::eval(const Vec_2<double> &pos, Vec_2<double> &v) const {
  const Vec_2<double>& u = get_disorder(pos);
#ifdef RANDOM_TORQUE
  v = Vec_2<double>(v.x * u.x - v.y * u.y, v.x * u.y + v.y * u.x);
#elif RANDOM_FIELD
  v += u;
#else
  if (v.dot(u) > 0) {
    v += u;
  } else {
    v -= u;
  }
#endif
}
