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

  template <class Par, class BoundayCondi>
  void interact(Par &p, const BoundayCondi &bc);
  template <class BoundaryCondi, class TRan>
  void move(double eta, double v0, TRan &myran, const BoundaryCondi &bc);
};

template <typename TRan>
VicsekPar_3::VicsekPar_3(TRan& myran, const Vec_3<double>& l,
                         const Vec_3<double>& origin) {
  const Vec_3<double> rand_vec3(myran.doub(), myran.doub(), myran.doub());
  pos = origin + rand_vec3 * l;
  sphere_point_picking(ori.x, ori.y, ori.z, myran);
  ori_next = ori;
}

template<class Par, class BoundayCondi>
void VicsekPar_3::interact(Par & p, const BoundayCondi & bc) {
  if(get_dis_square(this->pos, p.pos, bc) < 1) {
    ori_next += p.ori;
    ori += p.ori_next;
  }
}

template<class BoundaryCondi, class TRan>
void VicsekPar_3::move(double eta, double v0, TRan & myran,
                       const BoundaryCondi & bc) {
  ori_next.normalize();
  Vec_3<double> ax;
  ori_next.rand_perp_ax(ax, myran);
  ori_next.rotate((myran.doub() - 0.5) * eta, ax);
  ori = ori_next;
  pos += v0 * ori;
  bc.wrap(pos);
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