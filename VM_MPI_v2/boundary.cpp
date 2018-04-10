#include "boundary.h"

void BoundaryBase_2::ini(double Lx0, double Ly0, double x0, double y0) {
  Lx = Lx0;
  Ly = Ly0;
  x_min = x0;
  y_min = y0;
}

PBC_x_2::PBC_x_2(double Lx0, double Ly0, double x0, double y0) :
  BoundaryBase_2(Lx0, Ly0, x0, y0) {
  half_Lx = 0.5 * Lx;
  x_max = x_min + Lx;
}

void PBC_x_2::ini(double Lx0, double Ly0, double x0, double y0) {
  Lx = Lx0;
  Ly = Ly0;
  x_min = x0;
  y_min = y0;
  half_Lx = 0.5 * Lx;
  x_max = x_min + Lx;
}

PBC_xy_2::PBC_xy_2(double Lx0, double Ly0) :
  BoundaryBase_2(Lx0, Ly0) {
  half_Lx = 0.5 * Lx;
  half_Ly = 0.5 * Ly;
  x_max = x_min + Lx;
  y_max = y_min + Ly;
}

void PBC_xy_2::ini(double Lx0, double Ly0, double x0, double y0) {
  Lx = Lx0;
  Ly = Ly0;
  x_min = x0;
  y_min = y0;
  half_Lx = 0.5 * Lx;
  half_Ly = 0.5 * Ly;
  x_max = x_min + Lx;
  y_max = y_min + Ly;
}
