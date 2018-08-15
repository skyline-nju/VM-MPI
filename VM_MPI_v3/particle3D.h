#pragma once
#include "vect.h"
#include "domain3D.h"
class VicsekPar_3 {
public:
  VicsekPar_3() = default;
  VicsekPar_3(const Vec_3<double> &pos0, const Vec_3<double> &ori0)
    : pos(pos0), ori(ori0), ori_next(ori0) {}
  
  VicsekPar_3(const double* buf)
    : pos(buf[0], buf[1], buf[2]), ori(buf[3], buf[4], buf[5]), ori_next(ori) { std::cout << "bb" << std::endl; }
  template<typename TRan>
  VicsekPar_3(TRan &myran, const Vec_3<double> &l, const Vec_3<double> &origin);

  Vec_3<double> pos;
  Vec_3<double> ori;
  Vec_3<double> ori_next;

  template <class Par>
  void interact(Par &p);

  template <class Par>
  void interact(Par &p, const Domain_3 &domain);

  template <class Par>
  void interact(Par &p, const ExtDomain_3 &domain);

  template <class TRan>
  void move(double eta, double v0, TRan &myran);

  template <class TRan>
  void move(double eta, double v0, TRan &myran, const Domain_3 &domain);

  template <class TRan>
  void move(double eta, double v0, TRan &myran, const ExtDomain_3 &domain);
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
  }
}

template <class Par>
void VicsekPar_3::interact(Par& p, const Domain_3& domain) {
  Vec_3<double> dR = pos - p.pos;
  untangle_3(dR, domain.l(), domain.half_l());
  if (dR.square() < 1) {
    ori_next += p.ori;
    p.ori_next += ori;
  }
}

template <class Par>
void VicsekPar_3::interact(Par& p, const ExtDomain_3& domain) {
  Vec_3<double> dR = pos - p.pos;
  untangle_3(dR, domain.l(), domain.half_l(), domain.flag_ext());
  if (dR.square() < 1) {
    ori_next += p.ori;
    p.ori_next += ori;
  }
}

template <class TRan>
void VicsekPar_3::move(double eta, double v0, TRan& myran) {
  ori_next.normalize();
  Vec_3<double> ax{};
  ori_next.rand_perp_ax(ax, myran);
  ori_next.rotate((myran.doub() - 0.5) * eta * PI * 2, ax);
  ori = ori_next;
  pos += v0 * ori;
}

template <class TRan>
void VicsekPar_3::move(double eta, double v0, TRan& myran, const Domain_3& domain) {
  move(eta, v0, myran);
  tangle_3(pos, domain.origin(), domain.end_pnt(), domain.l());
}

template <class TRan>
void VicsekPar_3::move(double eta, double v0, TRan& myran, const ExtDomain_3& domain) {
  move(eta, v0, myran);
  tangle_3(pos, domain.origin(), domain.end_pnt(), domain.l(), domain.flag_ext());
}

template <class TPar>
void cal_order_para(const std::vector<TPar> &p_arr,
                    double &phi, Vec_3<double> &v_mean) {
  Vec_3<double> sum_v{};
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    sum_v += (*it).ori;
  }
  v_mean = sum_v / p_arr.size();
  phi = v_mean.module();
}

class VicsekPar_3_w_flag: public VicsekPar_3 {
public:
  VicsekPar_3_w_flag(): VicsekPar_3(), flag(false) {}
  VicsekPar_3_w_flag(const Vec_3<double> &pos, const Vec_3<double> &ori)
    : VicsekPar_3(pos, ori), flag(false) {}
  template<typename TRan>
  VicsekPar_3_w_flag(TRan &myran, const Vec_3<double> &l,
                     const Vec_3<double> &origin)
    : VicsekPar_3(myran, l, origin), flag(false) {}

  bool flag;
};

