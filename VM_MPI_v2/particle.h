#ifndef PARTICLE_H
#define PARTICLE_H
#include "rand.h"
#include "comn.h"
#include "vect.h"
#include "boundary.h"

/****************************************************************************/
/********************* Class for Vicske-style particle **********************/
/****************************************************************************/
class VicsekPar {
public:
  VicsekPar() {}
  VicsekPar(double x0, double y0, double vx0, double vy0) :
    x(x0), y(y0), vx(vx0), vy(vy0), vx_next(vx0), vy_next(vy0) {}
  VicsekPar(Ran &myran, double Lx, double Ly, double x0 = 0, double y0 = 0) {
    create_rand(myran, Lx, Ly, x0, y0);
  }
  
  void create_rand(Ran &myran, double Lx, double Ly, double x0, double y0);

  template <class Par, class BondaryCondi>
  void interact(Par &p, const BondaryCondi &bc);

  template <class Par, class BondaryCondi>
  void interact(Par &p, const BondaryCondi &bc, int &count);

  template <class BondaryCondi>
  void move(double eta, double v0, Ran &myran, const BondaryCondi &bc);

  template <class BondaryCondi>
  void move(double eta, double v0, Ran &myran, const BondaryCondi &bc,
    Vec_2<int> &d_cell);
  
  double x;
  double y;
  double vx;
  double vy;
  double vx_next;
  double vy_next;
};

inline void VicsekPar::create_rand(Ran &myran, double Lx, double Ly,
                              double x0, double y0) {
  x = x0 + myran.doub() * Lx;
  y = y0 + myran.doub() * Ly;
  double theta = myran.doub() * PI * 2;
  vx = vx_next = std::cos(theta);
  vy = vy_next = std::sin(theta);
}

template <class Par, class BondaryCondi>
inline void VicsekPar::interact(Par &p, const BondaryCondi &bc) {
  if (get_dis_square(*this, p, bc) < 1) {
    vx_next += p.vx;
    vy_next += p.vy;
    p.vx_next += vx;
    p.vy_next += vy;
  }
}

template <class Par, class BondaryCondi>
inline void VicsekPar::interact(Par &p, const BondaryCondi &bc, int &count) {
  if (get_dis_square(*this, p, bc) < 1) {
    vx_next += p.vx;
    vy_next += p.vy;
    p.vx_next += vx;
    p.vy_next += vy;
    count++;
  }
}

template <class BondaryCondi>
void VicsekPar::move(double eta, double v0, Ran &myran, const BondaryCondi &bc) {
  double tmp = std::sqrt(vx_next * vx_next + vy_next * vy_next);
  double c1 = vx_next / tmp;
  double s1 = vy_next / tmp;
  double noise = (myran.doub() - 0.5) * eta * 2 * PI;
  double c2 = std::cos(noise);
  double s2 = std::sin(noise);
  vx = vx_next = c1 * c2 - s1 * s2;
  vy = vy_next = c1 * s2 + c2 * s1;
  x += v0 * vx;
  y += v0 * vy;
  bc.wrap(*this);
}

template <class BondaryCondi>
void VicsekPar::move(double eta, double v0, Ran &myran, const BondaryCondi &bc,
  Vec_2<int> &d_cell) {
  double tmp = std::sqrt(vx_next * vx_next + vy_next * vy_next);
  double c1 = vx_next / tmp;
  double s1 = vy_next / tmp;
  double noise = (myran.doub() - 0.5) * eta * 2 * PI;
  double c2 = std::cos(noise);
  double s2 = std::sin(noise);
  vx = vx_next = c1 * c2 - s1 * s2;
  vy = vy_next = c1 * s2 + c2 * s1;
  d_cell.x = floor(x);
  d_cell.y = floor(y);
  x += v0 * vx;
  y += v0 * vy;
  d_cell.x = floor(x) - d_cell.x;
  d_cell.y = floor(y) - d_cell.y;
  bc.wrap(*this);
}

template <class _TPar>
void func_order_para(const std::vector<_TPar> &p_arr, double &phi, double &theta) {
  double svx = 0;
  double svy = 0;
  int nPar = p_arr.size();
  for (int i = 0; i < nPar; i++) {
    svx += p_arr[i].vx;
    svy += p_arr[i].vy;
  }
  phi = std::sqrt(svx * svx + svy * svy) / nPar;
  theta = std::atan2(svy, svx);
}
/* End of class VicsekPar */
#endif
