#include "cellList2D.h"

_CellListBase_2::_CellListBase_2(double Lx, double Ly, double x0, double y0,
  double r_cut, bool comm_x, bool comm_y) {
  if (comm_x == false && comm_y == false) {
    origin.x = 0;
    origin.y = 0;
    n_bins.x = int(Lx / r_cut);
    n_bins.y = int(Ly / r_cut);
    l_cell.x = Lx / n_bins.x;
    l_cell.y = Ly / n_bins.y;
    l_box.x = Lx;
    l_box.y = Ly;
  } else if (comm_x == false && comm_y == true) {
    n_bins.x = int(Lx / r_cut);
    n_bins.y = int(Ly / r_cut) + 2;
    l_cell.x = Lx / n_bins.x;
    l_cell.y = Ly / (n_bins.y - 2);
    l_box.x = Lx;
    l_box.y = Ly + 2 * l_cell.y;
    origin.x = x0;
    origin.y = y0 - l_cell.y;
  } else if (comm_x == true && comm_y == true) {
    n_bins.x = int(Lx / r_cut) + 2;
    n_bins.y = int(Ly / r_cut) + 2;
    l_cell.x = Lx / (n_bins.x - 2);
    l_cell.y = Ly / (n_bins.y - 2);
    l_box.x = Lx + 2 * l_cell.x;
    l_box.y = Ly + 2 * l_cell.y;
    origin.x = x0 - l_cell.x;
    origin.y = y0 - l_cell.y;
  }
  inverse_lc.x = 1 / l_cell.x;
  inverse_lc.y = 1 / l_cell.y;
  ncells = n_bins.x * n_bins.y;

  std::cout << "L_Box = " << l_box.x << "\t" << l_box.y << "\n"
    << "L_cell = " << l_cell.x << "\t" << l_cell.y << "\nx0 = "
    << origin.x << "\t y0 = " << origin.y << "\nnx = " << n_bins.x
    << "\tny = " << n_bins.y << "\n";
}
