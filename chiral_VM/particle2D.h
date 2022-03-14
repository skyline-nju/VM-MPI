#pragma once
#include "vect.h"
#include "domain2D.h"
#include "config.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

class Bird_2 {
public:
  Bird_2() = default;

  Bird_2(const Vec_2<double>& pos0, const Vec_2<double>& ori0)
    : pos(pos0), ori(ori0), tau(0.), n_neighb(1) {}

  Bird_2(const double* buf)
    : pos(buf[0], buf[1]), ori(buf[2], buf[3]), tau(0.), n_neighb(1) {}

  Bird_2(const float* buf)
    : pos(buf[0], buf[1]), ori(cos(buf[2]), sin(buf[2])), tau(0.), n_neighb(1) {}
  
  template <typename TRan>
  Bird_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& origin);

  void copy_from(const Vec_2<double>& pos_new, const Vec_2<double>& ori_new);

  void copy_to(double* dest, int& idx) const;

  double get_theta() const { return atan2(ori.y, ori.x); }

  Vec_2<double> pos;
  Vec_2<double> ori;   // phi = atan2(ori.y, ori.x)

  double tau;
  int n_neighb = 1;
};

template <typename TRan>
Bird_2::Bird_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& origin) {
  const Vec_2<double> rand_vec2(myran.doub(), myran.doub());
  pos = origin + rand_vec2 * l;
  circle_point_picking(ori.x, ori.y, myran);
  tau = 0.;
  n_neighb = 1;
}

inline void Bird_2::copy_from(const Vec_2<double>& pos_new,
                              const Vec_2<double>& ori_new) {
  pos = pos_new;
  ori = ori_new;
  tau = 0.;
  n_neighb = 1;
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
    double torque = p2.ori.y * p1.ori.x - p2.ori.x * p1.ori.y;
    p1.tau += torque;
    p2.tau -= torque;
    p1.n_neighb++;
    p2.n_neighb++;
  }
}

template <class Par, class TDomain>
void polar_align(Par& p1, Par& p2, const TDomain& domain) {
  Vec_2<double> dR = p2.pos - p1.pos;
  domain.untangle(dR);
  polar_align(p1, p2, dR);
}

//template <typename Par>
//void nematic_align(Par& p1, Par& p2, const Vec_2<double>& dR) {
//  // dR = p2.pos - p1.pos + offset
//  if (dR.square() < 1.) {
//    double sin_dtheta = p2.ori.y * p1.ori.x - p2.ori.x * p1.ori.y;
//    double cos_dtheta = p2.ori.x * p1.ori.x + p2.ori.y * p1.ori.y;
//    double torque = 2. * sin_dtheta * cos_dtheta;
//    p1.tau += torque;
//    p2.tau -= torque;
//
//    p1.n_neighb++;
//    p2.n_neighb++;
//  }
//}
//
//template <class Par, class TDomain>
//void nematic_align(Par& p1, Par& p2, const TDomain& domain) {
//  Vec_2<double> dR = p2.pos - p1.pos;
//  domain.untangle(dR);
//  nematic_align(p1, p2, dR);
//}

template <class Par>
void move_forward(Par& p, double v0, double dtheta) {
  const double c1 = p.ori.x;
  const double s1 = p.ori.y;
  //dtheta = scalar_noise + random_torque
  const double c2 = std::cos(dtheta);
  const double s2 = std::sin(dtheta);
  p.ori.x = c1 * c2 - s1 * s2;
  p.ori.y = c1 * s2 + c2 * s1;
  p.pos += v0 * p.ori;
}

template <class Par, class TDomain>
void move_forward(Par& p, double v0, double dtheta, const TDomain& domain) {
  move_forward(p, v0, dtheta);
  domain.tangle(p.pos);
}

