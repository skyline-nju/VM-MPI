#pragma once
#include "vect.h"
#include "domain3D.h"
#include "disorder.h"

class VicsekPar_3 {
public:
  VicsekPar_3() = default;

  VicsekPar_3(const Vec_3<double> &pos0, const Vec_3<double> &ori0)
    : pos(pos0), ori(ori0), ori_next(ori0) {}
  
  VicsekPar_3(const double* buf)
    : pos(buf[0], buf[1], buf[2]), ori(buf[3], buf[4], buf[5]), ori_next(ori) {}

  template<typename TRan>
  VicsekPar_3(TRan &myran, const Vec_3<double> &l, const Vec_3<double> &origin);


  template <class Par>
  void interact(Par &p);

  template <class Par>
  void interact(Par &p, const Vec_3<double> &offset);

  template <class Par>
  void interact(Par &p, const Domain_3 &domain);

  template<class Par>
  void interact_w_ghost(Par &ghost, const Domain_3 &domain);

  template <class TRan>
  void move(double eta, double v0, TRan &myran);

  template <class TRan>
  void move(double eta, double v0, const RandTorque &torque, TRan &myran);

  template <class TRan>
  void move(double eta, double v0, TRan &myran, const Domain_3 &domain);

  template <class TRan>
  void move(double eta, double v0, const RandTorque &torque,  TRan &myran, const Domain_3 &domain);

  void copy(double* dest, int& idx) const;

  Vec_3<double> pos;
  Vec_3<double> ori;
  Vec_3<double> ori_next;
#ifdef RAND_FIELD
  int n_neighb=1;
#endif
};

template <typename TRan>
VicsekPar_3::VicsekPar_3(TRan& myran, const Vec_3<double>& l,
                         const Vec_3<double>& origin): pos(), ori(), ori_next() {
  const Vec_3<double> rand_vec3(myran.doub(), myran.doub(), myran.doub());
  pos = origin + rand_vec3 * l;
  sphere_point_picking(ori.x, ori.y, ori.z, myran);
  ori_next = ori;
}

template <class Par>
void VicsekPar_3::interact(Par& p) {
  Vec_3<double> dR = pos - p.pos;
  if (dR.square() < 1) {
    ori_next += p.ori;
    p.ori_next += ori;
#ifdef RAND_FIELD
  n_neighb += 1;
  p.n_neibhb += 1;
#endif
  }
}

template <class Par>
void VicsekPar_3::interact(Par& p, const Vec_3<double>& offset) {
  Vec_3<double> dR = p.pos - pos + offset;
  if (dR.square() < 1) {
    ori_next += p.ori;
    p.ori_next += ori;
#ifdef RAND_FIELD
  n_neighb += 1;
  p.n_neibhb += 1;
#endif
  }
}


template <class Par>
void VicsekPar_3::interact(Par& p, const Domain_3& domain) {
  Vec_3<double> dR = pos - p.pos;
  untangle_3(dR, domain.gl_l(), domain.gl_half_l(), domain.flag_comm());
  if (dR.square() < 1) {
    ori_next += p.ori;
    p.ori_next += ori;
#ifdef RAND_FIELD
  n_neighb += 1;
  p.n_neibhb += 1;
#endif
  }
}

template<class Par>
void VicsekPar_3::interact_w_ghost(Par & ghost, const Domain_3 & domain) {
  Vec_3<double> dR = pos - ghost.pos;
  untangle_3(dR, domain.gl_l(), domain.gl_half_l(), domain.flag_comm());
  if (dR.square() < 1) {
    ori_next += ghost.ori;
    ghost.ori_next += ori;
#ifdef RAND_FIELD
  n_neighb += 1;
  p.n_neibhb += 1;
#endif
  }
}

template <class TRan>
void VicsekPar_3::move(double eta, double v0, TRan& myran) {
  ori_next.normalize();
  double theta = (myran.doub() - 0.5) * eta * PI * 2;
  ori_next.rotate_rand(theta, myran);
  ori = ori_next;
  pos += v0 * ori;
#ifdef RAND_FIELD
  n_neighb = 1;
#endif
}

template <class TRan>
void VicsekPar_3::move(double eta, double v0, const RandTorque& torque, TRan& myran) {
  ori_next.normalize();
  double theta = (myran.doub() - 0.5) * eta * PI * 2;
  ori_next.rotate_rand(theta, myran);
  torque.rotate(ori_next);
  ori = ori_next;
  pos += v0 * ori;
}

template <class TRan>
void VicsekPar_3::move(double eta, double v0, TRan& myran, const Domain_3& domain) {
  move(eta, v0, myran);
  tangle_3(pos, domain.gl_l(), domain.flag_comm());
}

template <class TRan>
void VicsekPar_3::move(double eta, double v0, const RandTorque& torque, TRan& myran, const Domain_3& domain) {
  move(eta, v0, torque, myran);
  tangle_3(pos, domain.gl_l(), domain.flag_comm());
}


inline void VicsekPar_3::copy(double *dest, int &idx) const {
  dest[idx]     = pos.x;
  dest[idx + 1] = pos.y;
  dest[idx + 2] = pos.z;
  dest[idx + 3] = ori.x;
  dest[idx + 4] = ori.y;
  dest[idx + 5] = ori.z;
  idx += 6;
}

