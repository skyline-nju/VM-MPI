#include "boundary3D.h"
#include "comn.h"

PBC_xy_3::PBC_xy_3(const Vec3d& l, const Vec3d& origin)
  : BoundaryBase_3(l, origin) {
  half_lx_ = l.x * 0.5;
  half_ly_ = l.y * 0.5;
  x_max_ = origin.x + l.x;
  y_max_ = origin.y + l.y;
}

void PBC_xy_3::ini(const Vec3d & l, const Vec3d & origin) {
  l_ = l;
  origin_ = origin;
  half_lx_ = l.x * 0.5;
  half_ly_ = l.y * 0.5;
  x_max_ = origin.x + l.x;
  y_max_ = origin.y + l.y;
}

void PBC_xy_3::nearest_dis(Vec3d& dis) const {
  untangle_1(dis.x, l_.x, half_lx_);
  untangle_1(dis.y, l_.y, half_ly_);
}

PBC_xyz_3::PBC_xyz_3(const Vec3d& l, const Vec3d& origin)
  : BoundaryBase_3(l, origin) {
  half_l_ = l_ * 0.5;
  end_ = origin_ + l_;
}

void PBC_xyz_3::ini(const Vec3d& l, const Vec3d& origin) {
  l_ = l;
  origin_ = origin;
  half_l_ = l_ * 0.5;
  end_ = origin_ + l_;
}

void PBC_xyz_3::nearest_dis(Vec3d& dis) const {
  untangle_1(dis.x, l_.x, half_l_.x);
  untangle_1(dis.y, l_.y, half_l_.y);
  untangle_1(dis.z, l_.z, half_l_.z);
}
