#pragma once
#include "vect.h"
#include "domain2D.h"
#include "config.h"
#include "disorder2D.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

class Bird_2 {
public:
  Bird_2() = default;
  Bird_2(const Vec_2<double>& pos0, const Vec_2<double>& ori0) : pos(pos0), ori(ori0), ori_next(ori0) {}
  Bird_2(const double* buf) : pos(buf[0], buf[1]), ori(buf[2], buf[3]), ori_next(ori) {}
  template <typename TRan>
  Bird_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& origin);

  void copy_from(const Vec_2<double>& pos_new, const Vec_2<double>& ori_new);

  void copy_to(double* dest, int& idx) const;

  double theta() const { return atan2(ori.y, ori.x); }

  Vec_2<double> pos;
  Vec_2<double> ori;
  Vec_2<double> ori_next;

#if defined RANDOM_FIELD || defined RANDOM_OBSTACLE
  int n_neighb = 1;
#endif
};

template <typename TRan>
Bird_2::Bird_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& origin) {
  const Vec_2<double> rand_vec2(myran.doub(), myran.doub());
  pos = origin + rand_vec2 * l;
  circle_point_picking(ori.x, ori.y, myran);
  ori_next = ori;
}

inline void Bird_2::copy_from(const Vec_2<double>& pos_new, const Vec_2<double>& ori_new) {
  pos = pos_new;
  ori = ori_next = ori_new;
#if defined RANDOM_FIELD || defined RANDOM_OBSTACLE
  n_neighb = 1;
#endif
}

inline void Bird_2::copy_to(double* dest, int& idx) const {
  dest[idx] = pos.x;
  dest[idx + 1] = pos.y;
  dest[idx + 2] = ori.x;
  dest[idx + 3] = ori.y;
  idx += 4;
}

template <typename Par>
void polar_align(Par &p1, Par &p2, const Vec_2<double>& dR) {
  // dR = p2.pos - p1.pos + offset
  if (dR.square() < 1.) {
    p1.ori_next += p2.ori;
    p2.ori_next += p1.ori;
#if defined RANDOM_FIELD || defined RANDOM_OBSTACLE
    p1.n_neighb++;
    p2.n_neighb++;
#endif
  }
}

template <class Par, class TDomain>
void polar_align(Par& p1, Par& p2, const TDomain& domain) {
  Vec_2<double> dR = p2.pos - p1.pos;
  domain.untangle(dR);
  polar_align(p1, p2, dR);
}

template <typename Par>
void nematic_align(Par& p1, Par& p2, const Vec_2<double>& dR) {
  // dR = p2.pos - p1.pos + offset
  if (dR.square() < 1.) {
    if (p1.ori.dot(p2.ori) > 0) {
      p1.ori_next += p2.ori;
      p2.ori_next += p1.ori;
    } else {
      p1.ori_next -= p2.ori;
      p2.ori_next -= p1.ori;
    }
#if defined RANDOM_FIELD || defined RANDOM_OBSTACLE
    p1.n_neighb++;
    p2.n_neighb++;
#endif
  }
}

template <class Par, class TDomain>
void nematic_align(Par& p1, Par& p2, const TDomain& domain) {
  Vec_2<double> dR = p2.pos - p1.pos;
  domain.untangle(dR);
  nematic_align(p1, p2, dR);
}

template <class Par>
void move_forward(Par& p, double v0, double dtheta) {
  p.ori_next.normalize();
  const double c1 = p.ori_next.x;
  const double s1 = p.ori_next.y;

  //dtheta = scalar_noise + random_torque
  const double c2 = std::cos(dtheta);
  const double s2 = std::sin(dtheta);
  p.ori.x = p.ori_next.x = c1 * c2 - s1 * s2;
  p.ori.y = p.ori_next.y = c1 * s2 + c2 * s1;
  p.pos += v0 * p.ori;
#if defined RANDOM_FIELD || defined RANDOM_OBSTACLE
  p.n_neighb = 1;
#endif
}

template <class Par, class TDomain>
void move_forward(Par& p, double v0, double dtheta, const TDomain& domain) {
  move_forward(p, v0, dtheta);
  domain.tangle(p.pos);
}


template <typename TPar>
void get_vel_sum(double* vel_sum, const std::vector<TPar>& p_arr) {
  vel_sum[0] = vel_sum[1] = 0.;
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
#ifdef POLAR_ALIGN
    vel_sum[0] += (*it).ori.x;
    vel_sum[1] += (*it).ori.y;
#else
    vel_sum[0] += (*it).ori.x * (*it).ori.x - (*it).ori.y * (*it).ori.y;
    vel_sum[1] += 2 * (*it).ori.x * (*it).ori.y;
#endif
  }
}

template <typename TPar>
void get_vel_mean(double* vel_mean, const std::vector<TPar>& p_arr) {
  get_vel_sum(vel_mean, p_arr);
  vel_mean[0] /= p_arr.size();
  vel_mean[1] /= p_arr.size();
}

#ifdef USE_MPI
template <typename TPar>
void get_mean_vel(double* vel_mean, const std::vector<TPar>& p_arr,
                  int gl_np, bool flag_broadcast, MPI_Comm group_comm) {
  get_vel_sum(vel_mean, p_arr);
  double gl_vel_sum[2];
  MPI_Reduce(vel_mean, gl_vel_sum, 2, MPI_DOUBLE, MPI_SUM, 0, group_comm);
  int my_rank;
  MPI_Comm_rank(group_comm, &my_rank);
  if (my_rank == 0) {
    vel_mean[0] = gl_vel_sum[0] / gl_np;
    vel_mean[1] = gl_vel_sum[1] / gl_np;
  }
  if (flag_broadcast) {
    MPI_Bcast(vel_mean, 2, MPI_DOUBLE, 0, group_comm);
  }
}
#endif