#pragma once
#include "vect.h"
#include "domain2D.h"
#include "config.h"
#include "disorder2D.h"

class Bird_2 {
public:
  Bird_2() = default;
  Bird_2(const Vec_2<double>& pos0, const Vec_2<double>& ori0) : pos(pos0), ori(ori0), ori_next(ori0) {}
  Bird_2(const double* buf) : pos(buf[0], buf[1]), ori(buf[2], buf[3]), ori_next(ori) {}
  template <typename TRan>
  Bird_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& origin);


  void copy_from(const Vec_2<double>& pos_new, const Vec_2<double>& ori_new);

  void copy_to(double* dest, int& idx) const;

  Vec_2<double> pos;
  Vec_2<double> ori;
  Vec_2<double> ori_next;

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
  }
}

template <class Par>
void polar_align(Par& p1, Par& p2, const Domain_2& domain) {
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
  }
}

template <class Par>
void nematic_align(Par& p1, Par& p2, const Domain_2& domain) {
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
}

template <class Par>
void move_forward(Par& p, double v0, double dtheta, const Domain_2& domain) {
  move_forward(p, v0, dtheta);
  domain.tangle(p.pos);
}


/*

class VicsekPar_2 {
public:
  VicsekPar_2() = default;

  VicsekPar_2(const Vec_2<double> &pos0, const Vec_2<double> &ori0)
    : pos(pos0), ori(ori0), ori_next(ori0) {}

  VicsekPar_2(const double* buf)
    : pos(buf[0], buf[1]), ori(buf[2], buf[3]), ori_next(ori) {}

  template<typename TRan>
  VicsekPar_2(TRan &myran, const Vec_2<double> &l, const Vec_2<double> &origin);

  void copy_from(const Vec_2<double> &pos_new, const Vec_2<double> &ori_new);

  void copy_to(double* dest, int& idx) const;

  template <class Par>
  void interact(Par &p);

  template <class Par>
  void interact(const Vec_2<double> &dR, Par &p);

  template <class Par>
  void interact(Par &p, const Vec_2<double> &offset) { interact(p.pos - pos + offset, p); }

  template <class Par>
  void interact(Par &p, const Domain_2 &domain);

  template <class TRan>
  void move(double eta, double v0, TRan &myran);

  template <class TRan>
  void move(double eta, double v0, TRan &myran, const Domain_2 &domain);

  template <class TRan>
  void move(double eta, double v0, TRan &myran, const Domain_2 &domain, const Disorder_2 &disorder);

#ifdef DENSITY_NOISE
  template <class TRan>
  void move_density_noise(double eta, double v0, double eps, TRan &myran, const Domain_2 &domain);
#endif

  Vec_2<double> pos;
  Vec_2<double> ori;
  Vec_2<double> ori_next;
#ifdef COUNT_NEIGHBOR
  int n_neighbor = 0;
#endif
};

template <typename TRan>
VicsekPar_2::VicsekPar_2(TRan& myran, const Vec_2<double>& l,
                         const Vec_2<double>& origin) : pos(), ori(), ori_next() {
  const Vec_2<double> rand_vec2(myran.doub(), myran.doub());
  pos = origin + rand_vec2 * l;
  circle_point_picking(ori.x, ori.y, myran);
  ori_next = ori;
}

inline void VicsekPar_2::copy_from(const Vec_2<double>& pos_new,
  const Vec_2<double>& ori_new) {
  pos = pos_new;
  ori = ori_next = ori_new;
}

inline void VicsekPar_2::copy_to(double* dest, int& idx) const {
  dest[idx] = pos.x;
  dest[idx + 1] = pos.y;
  dest[idx + 2] = ori.x;
  dest[idx + 3] = ori.y;
  idx += 4;
}

template <class Par>
void VicsekPar_2::interact(Par& p) {
  interact(pos - p.pos, p);
}

template <class Par>
void VicsekPar_2::interact(const Vec_2<double> &dR, Par &p) {
  if (dR.square() < 1.) {
#ifdef POLAR_ALIGN
    ori_next += p.ori;
    p.ori_next += ori;
#else
// nematic aligning
    if (ori.dot(p.ori) > 0) {
      ori_next += p.ori;
      p.ori_next += ori;
    } else {
      ori_next -= p.ori;
      p.ori_next -= ori;
    }
#endif
#ifdef COUNT_NEIGHBOR
    n_neighbor++;
    p.n_neighbor += 1;
#endif
  }
}

template <class Par>
void VicsekPar_2::interact(Par& p, const Domain_2& domain) {
  Vec_2<double> dR = pos - p.pos;
  if (!domain.flag_comm().x) {
    untangle_1(dR.x, domain.gl_l().x, domain.gl_half_l().x);
  }
  if (!domain.flag_comm().y) {
    untangle_1(dR.y, domain.gl_l().y, domain.gl_half_l().y);
  }
  interact(dR, p);
}

template <class TRan>
void VicsekPar_2::move(double eta, double v0, TRan& myran) {
#ifdef SCALAR_NOISE
  ori_next.normalize();
  const double c1 = ori_next.x;
  const double s1 = ori_next.y;
 
  const double theta = (myran.doub() - 0.5) * eta * PI * 2;
  const double c2 = std::cos(theta);
  const double s2 = std::sin(theta);
  ori.x = ori_next.x = c1 * c2 - s1 * s2;
  ori.y = ori_next.y = c1 * s2 + c2 * s1;
#else
  const double ori_next_module = ori_next.module();
  Vec_2<double> rand_unit_vec{};
  circle_point_picking(rand_unit_vec.x, rand_unit_vec.y, myran);
  ori_next += eta * n_neighbor * rand_unit_vec;
  ori_next.normalize();
  ori = ori_next;
  n_neighbor = 1;
#endif
  pos += v0 * ori;
}

template <class TRan>
void VicsekPar_2::move(double eta, double v0, TRan& myran, const Domain_2& domain) {
  move(eta, v0, myran);
  domain.tangle(pos);
}

template <class TRan>
void VicsekPar_2::move(double eta, double v0, TRan& myran,
                       const Domain_2& domain, const Disorder_2& disorder) {
#ifdef SCALAR_NOISE
#ifdef DISORDER_ON
#ifndef RANDOM_TORQUE
  ori_next /= (n_neighbor + 1.0);
  disorder.eval(pos, ori_next);
#endif
#endif
  ori_next.normalize();
  const double c1 = ori_next.x;
  const double s1 = ori_next.y;
  const double theta = (myran.doub() - 0.5) * eta * PI * 2;
  const double c2 = std::cos(theta);
  const double s2 = std::sin(theta);
  ori_next.x = c1 * c2 - s1 * s2;
  ori_next.y = c1 * s2 + c2 * s1;
#ifdef RANDOM_TORQUE
  disorder.eval(pos, ori_next);
#endif
#else
  Vec_2<double> rand_unit_vec{};
  circle_point_picking(rand_unit_vec.x, rand_unit_vec.y, myran);
  ori_next += eta * (n_neighbor + 1) * rand_unit_vec;
#ifndef RANDOM_TORQUE
  disorder.eval(pos, ori_next, n_neighbor + 1);
#endif
  ori_next.normalize();
#ifdef RANDOM_TORQUE
  disorder.eval(pos, ori_next);
#endif
#endif
#ifdef COUNT_NEIGHBOR
  n_neighbor = 0;
#endif
  ori = ori_next;
  pos += v0 * ori;
  domain.tangle(pos);
}

#ifdef DENSITY_NOISE
template <class TRan>
void VicsekPar_2::move_density_noise(double eta, double v0, double eps, TRan &myran, const Domain_2& domain) {
#ifdef SCALAR_NOISE
  ori_next.normalize();
  const double c1 = ori_next.x;
  const double s1 = ori_next.y;

  const double theta = (myran.doub() - 0.5) * eta * PI * 2 + n_neighbor * (myran.doub() - 0.5) * eps * PI * 2;
  const double c2 = std::cos(theta);
  const double s2 = std::sin(theta);
  ori.x = ori_next.x = c1 * c2 - s1 * s2;
  ori.y = ori_next.y = c1 * s2 + c2 * s1;
  n_neighbor = 0;
#else
  //todo
  Vec_2<double> rand_unit_vec{};
  circle_point_picking(rand_unit_vec.x, rand_unit_vec.y, myran);
  ori_next += eta * n_neighbor * rand_unit_vec;
  ori_next.normalize();
  ori_next.rotate(torque);
  ori = ori_next;
  n_neighbor = 0;
#endif
  pos += v0 * ori;
  domain.tangle(pos);
}
#endif
*/