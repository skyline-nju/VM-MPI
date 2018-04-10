#ifndef BOUNDARY_H
#define BOUNDARY_H
#include "vect.h"
#include <iostream>

template <typename _T1, typename _T2, typename _TBC>
inline double get_dis_square(const _T1 &p1, const _T2 &p2, const _TBC &bc) {
  Vec_2<double> dis(p1.x - p2.x, p1.y - p2.y);
  bc.nearest_dis(dis);
  return dis.square();
}

class BoundaryBase_2 {
public:
  BoundaryBase_2() {}
  BoundaryBase_2(double Lx0, double Ly0, double x0 = 0, double y0 = 0) :
    Lx(Lx0), Ly(Ly0), x_min(x0), y_min(y0) {}
  void ini(double Lx0, double Ly0, double x0=0, double y0=0);
  double get_Lx() const { return Lx; }
  double get_Ly() const { return Ly; }
protected:
  double Lx;
  double Ly;
  double x_min;
  double y_min;
};

class PBC_x_2 : public BoundaryBase_2 {
public:
  PBC_x_2() : BoundaryBase_2() {}
  PBC_x_2(double Lx0, double Ly0, double x0=0, double y0=0);
  void ini(double Lx0, double Ly0, double x0 = 0, double y0 = 0);
  void nearest_dis(Vec_2<double> &dis) const;
  template <typename T>
  void wrap(T &p) const;

protected:
  double half_Lx;
  double x_max;

};

inline void PBC_x_2::nearest_dis(Vec_2<double> &dis) const {
  if (dis.x < -half_Lx) {
    dis.x += Lx;
  } else if (dis.x > half_Lx) {
    dis.x -= Lx;
  }
}

template <typename T>
void PBC_x_2::wrap(T &p) const {
  if (p.x < x_min) {
    p.x += Lx;
  } else if (p.x >= x_max) {
    p.x -= Lx;
  }
}

class PBC_xy_2 : public BoundaryBase_2 {
public:
  PBC_xy_2() : BoundaryBase_2() {}
  PBC_xy_2(double Lx0, double Ly0);
  void ini(double Lx0, double Ly0, double x0 = 0, double y0 = 0);
  void nearest_dis(Vec_2<double> &dis) const;
  template <typename T>
  void wrap(T &p) const;

protected:
  double half_Lx;
  double half_Ly;
  double x_max;
  double y_max;
};

inline void PBC_xy_2::nearest_dis(Vec_2<double> &dis) const {
  if (dis.x < -half_Lx) {
    dis.x += Lx;
  } else if (dis.x > half_Lx) {
    dis.x -= Lx;
  }
  if (dis.y < -half_Ly) {
    dis.y += Ly;
  } else if (dis.y > half_Ly) {
    dis.y -= Ly;
  }
}

template<typename T>
inline void PBC_xy_2::wrap(T & p) const {
  if (p.x < x_min) {
    p.x += Lx;
  } else if (p.x >= x_max) {
    p.x -= Lx;
  }
  if (p.y < x_min) {
    p.y += Ly;
  } else if (p.y >= x_max) {
    p.y -= Ly;
  }
}
#endif

