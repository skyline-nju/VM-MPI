#pragma once
#include "vect.h"
#include "domain2D.h"
#include "config.h"

class VicsekPar_2 {
public:
  VicsekPar_2() = default;

  VicsekPar_2(const Vec_2<double> &pos0, const Vec_2<double> &ori0)
    : pos(pos0), ori(ori0), ori_next(ori0) {}

  VicsekPar_2(const double* buf)
    : pos(buf[0], buf[1]), ori(buf[2], buf[3]), ori_next(ori) {}

  template<typename TRan>
  VicsekPar_2(TRan &myran, const Vec_2<double> &l, const Vec_2<double> &origin);

  template <class Par>
  void interact(Par &p);

  template <class Par>
  void interact(Par &p, const Vec_2<double> &offset);

  template <class Par>
  void interact(Par &p, const Domain_2 &domain);

  template <class TRan>
  void move(double eta, double v0, TRan &myran);

  template <class TRan>
  void move(double eta, double v0, TRan &myran, const Domain_2 &domain);

  template <class TRan>
  void move_torque(double eta, double v0, double torque, TRan &myran);

  template <class TRan>
  void move_torque(double eta, double v0, double torque, TRan &myran, const Domain_2 &domain);

  template <class TRan>
  void move_field(double eta, double v0, double h, TRan &myran);

  template <class TRan>
  void move_field(double eta, double v0, double h, TRan &myran, const Domain_2 &domain);

  template <class TRan>
  void move_density_noise(double eta, double v0, double eps, TRan &myran, const Domain_2 &domain);

  void copy(double* dest, int& idx) const;

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

template <class Par>
void VicsekPar_2::interact(Par& p) {
  Vec_2<double> dR = pos - p.pos;
  if (dR.square() < 1) {
    ori_next += p.ori;
    p.ori_next += ori;
#ifdef COUNT_NEIGHBOR
    n_neighbor++;
    p.n_neighbor++;
#endif
  }
}

template <class Par>
void VicsekPar_2::interact(Par& p, const Vec_2<double>& offset) {
  Vec_2<double> dR = p.pos - pos + offset;
  if (dR.square() < 1) {
    ori_next += p.ori;
    p.ori_next += ori;
#ifdef COUNT_NEIGHBOR
    n_neighbor++;
    p.n_neighbor++;
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

  if (dR.square() < 1) {
    ori_next += p.ori;
    p.ori_next += ori;
#ifdef COUNT_NEIGHBOR
    n_neighbor++;
    p.n_neighbor++;
#endif
  }
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
void VicsekPar_2::move_torque(double eta, double v0, double torque, TRan& myran) {
#ifdef SCALAR_NOISE
  ori_next.normalize();
  const double c1 = ori_next.x;
  const double s1 = ori_next.y;

  const double theta = (myran.doub() - 0.5) * eta * PI * 2 + torque;
  const double c2 = std::cos(theta);
  const double s2 = std::sin(theta);
  ori.x = ori_next.x = c1 * c2 - s1 * s2;
  ori.y = ori_next.y = c1 * s2 + c2 * s1;
#else
  Vec_2<double> rand_unit_vec{};
  circle_point_picking(rand_unit_vec.x, rand_unit_vec.y, myran);
  ori_next += eta * n_neighbor * rand_unit_vec;
  ori_next.normalize();
  ori_next.rotate(torque);
  ori = ori_next;
  n_neighbor = 1;
#endif
  pos += v0 * ori;
}

template <class TRan>
void VicsekPar_2::move_torque(double eta, double v0, double torque, TRan& myran, const Domain_2& domain) {
  move_torque(eta, v0, torque, myran);
  domain.tangle(pos);
}

template <class TRan>
void VicsekPar_2::move_field(double eta, double v0, double h, TRan& myran) {
  const double ori_next_module = ori_next.module();
  ori_next.x += h * ori_next_module;
#ifdef SCALAR_NOISE
  // todo
  // ori_next.normalize();
  // const double c1 = ori_next.x;
  // const double s1 = ori_next.y;
  // const double theta = (myran.doub() - 0.5) * eta * PI * 2 + torque;
  // const double c2 = std::cos(theta);
  // const double s2 = std::sin(theta);
  // ori.x = ori_next.x = c1 * c2 - s1 * s2;
  // ori.y = ori_next.y = c1 * s2 + c2 * s1;
#else
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
void VicsekPar_2::move_field(double eta, double v0, double h, TRan& myran, const Domain_2& domain) {
  move_field(eta, v0, h, myran);
  domain.tangle(pos);
}

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

inline void VicsekPar_2::copy(double *dest, int &idx) const {
  dest[idx] = pos.x;
  dest[idx + 1] = pos.y;
  dest[idx + 2] = ori.x;
  dest[idx + 3] = ori.y;
  idx += 4;
}
