#pragma once
#include "vect.h"
#include "boundary3D.h"
class VicsekPar_3 {
public:
  VicsekPar_3() = default;
  VicsekPar_3(const Vec_3<double> &pos0, const Vec_3<double> &ori0)
    : pos(pos0), ori(ori0), ori_next(ori0) {}
  template<typename TRan>
  VicsekPar_3(TRan &myran, const Vec_3<double> &l, const Vec_3<double> &origin);

  Vec_3<double> pos;
  Vec_3<double> ori;
  Vec_3<double> ori_next;

  template <class Par>
  void interact(Par &p);

  template <class Par, class BoundayCondi>
  void interact(Par &p, const BoundayCondi &bc);

  template <class Par>
  void interact(Par &p, const Vec_3<double> &l, const Vec_3<double> &half_l);


  template <class BoundaryCondi, class TRan>
  void my_move(double eta, double v0, TRan &myran, const BoundaryCondi &bc);

  template <class TRan>
  void my_move(double eta, double v0, TRan &myran, const Vec_3<double> &origin,
    const Vec_3<double> &end, const Vec_3<double> &l);

  template <class BoundaryCondi, class TRan>
  void move(double eta, double v0, TRan &myran, const BoundaryCondi &bc) {
    my_move(eta, v0, myran, bc);
  }

  template <class TRan>
  void move(double eta, double v0, TRan &myran, const Vec_3<double> &origin,
    const Vec_3<double> &end, const Vec_3<double> &l) {
    my_move(eta, v0, myran, origin, end, l);
  }
};

template <typename TRan>
VicsekPar_3::VicsekPar_3(TRan& myran, const Vec_3<double>& l,
                         const Vec_3<double>& origin) {
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

template<class Par, class BoundayCondi>
void VicsekPar_3::interact(Par & p, const BoundayCondi & bc) {
  if(get_dis_square(this->pos, p.pos, bc) < 1) {
    ori_next += p.ori;
    p.ori_next += ori;
  }
}

template<class Par>
void VicsekPar_3::interact(Par & p, const Vec_3<double>& l, const Vec_3<double>& half_l) {
  Vec_3<double> dR = pos - p.pos;
  untangle_1(dR.x, l.x, half_l.x);
  untangle_1(dR.y, l.y, half_l.y);
  untangle_1(dR.z, l.z, half_l.z);
  if (dR.square() < 1) {
    ori_next += p.ori;
    p.ori_next += ori;
  }
}

template<class BoundaryCondi, class TRan>
void VicsekPar_3::my_move(double eta, double v0, TRan & myran,
                          const BoundaryCondi & bc) {
  ori_next.normalize();
  Vec_3<double> ax;
  ori_next.rand_perp_ax(ax, myran);
  ori_next.rotate((myran.doub() - 0.5) * eta * PI * 2, ax);
  ori = ori_next;
  pos += v0 * ori;
  bc.wrap(pos);
}

template <class TRan>
void VicsekPar_3::my_move(double eta, double v0, TRan& myran, const Vec_3<double>& origin, const Vec_3<double>& end,
  const Vec_3<double>& l) {
  ori_next.normalize();
  Vec_3<double> ax;
  ori_next.rand_perp_ax(ax, myran);
  ori_next.rotate((myran.doub() - 0.5) * eta * PI * 2, ax);
  ori = ori_next;
  pos += v0 * ori;
  tangle_1(pos.x, origin.x, end.x, l.x);
  tangle_1(pos.y, origin.y, end.y, l.y);
  tangle_1(pos.z, origin.z, end.z, l.z);
}

template <class TPar>
void cal_order_para(const std::vector<TPar> &p_arr,
                    double &phi, Vec_3<double> &v_mean) {
  Vec_3<double> sum_v;
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
  VicsekPar_3_w_flag(TRan &myran, const Vec_3<double> &l, const Vec_3<double> &origin)
    : VicsekPar_3(myran, l, origin), flag(false) {}

  template <class TRan, class BoundaryCondi>
  void move(double eta, double v0, TRan &myran, const BoundaryCondi &bc);
  
  template <class TRan>
  void move(double eta, double v0, TRan &myran, const Vec_3<double> &origin,
    const Vec_3<double> &end, const Vec_3<double> &l);

  bool flag;
};

template <class TRan, class BoundaryCondi>
void VicsekPar_3_w_flag::move(double eta, double v0, TRan& myran,
                              const BoundaryCondi& bc) {
  my_move(eta, v0, myran, bc);
  flag = true;
}

template <class TRan>
void VicsekPar_3_w_flag::move(double eta, double v0, TRan& myran, const Vec_3<double>& origin, const Vec_3<double>& end,
  const Vec_3<double>& l) {
  my_move(eta, v0, myran, origin, end, l);
  flag = true;
}
