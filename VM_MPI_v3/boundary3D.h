#pragma once

#include "vect.h"
#include <comn.h>

template <typename T1, typename T2, typename BoundaryCondi>
double get_dis_square(const T1& p1, const T2& p2, const BoundaryCondi &bc) {
  Vec_3<double> dis(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
  bc.nearest_dis(dis);
  return dis.square();
}

class BoundaryBase_3 {
public:
  typedef Vec_3<double> Vec3d;
  BoundaryBase_3() = default;
  explicit BoundaryBase_3(const Vec3d &l, const Vec3d &origin=Vec3d())
    : l_(l), origin_(origin){}
  double get_Lx() const { return l_.x; }
  double get_Ly() const { return l_.y; }
  double get_Lz() const { return l_.z; }

protected:
  Vec_3<double> l_;
  Vec_3<double> origin_;
};

class PBC_xy_3: public BoundaryBase_3 {
public:
  PBC_xy_3() = default;
  explicit PBC_xy_3(const Vec3d &l, const Vec3d &origin = Vec3d());
  void ini(const Vec3d &l, const Vec3d &origin = Vec3d());
  void nearest_dis(Vec3d &dis) const;
  template <typename T>
  void wrap(T &p) const;

protected:
  double half_lx_;
  double half_ly_;
  double x_max_;
  double y_max_;
};

template <typename T>
void PBC_xy_3::wrap(T& p) const {
  tangle_1(p.x, origin_.x, x_max_, l_.x);
  tangle_1(p.y, origin_.y, y_max_, l_.y);
}

class PBC_xyz_3: public BoundaryBase_3 {
public:
  PBC_xyz_3() = default;
  explicit PBC_xyz_3(const Vec3d &l, const Vec3d &origin = Vec3d());
  void ini(const Vec3d &l, const Vec3d &origin = Vec3d());
  void nearest_dis(Vec3d &dis) const;
  template <typename T>
  void wrap(T &p) const;

protected:
  Vec3d half_l_;
  Vec3d end_;
};

template <typename T>
void PBC_xyz_3::wrap(T& p) const {
  tangle_1(p.x, origin_.x, end_.x, l_.x);
  tangle_1(p.y, origin_.y, end_.y, l_.y);
  tangle_1(p.z, origin_.z, end_.z, l_.z);
}
