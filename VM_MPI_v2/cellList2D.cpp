#include "cellList2D.h"
#include <set>

const int CellListBase_2::cell_offset[4][2] = {
  // {offset_y, offset_x}
  { 0,  1 },
  { 1, -1 },
  { 1,  0 },
  { 1,  1 }
};


CellListBase_2::CellListBase_2(const Vec2d & l,
                               double r_cut,
                               const Vec2d &gl_l,
                               const Vec2d & origin,
                               const Vec_2<bool>& flag_ext)
                               : n_(), origin_(origin), lc_recip_(), l_(l), gl_l_(gl_l),
                               flag_ext_(flag_ext) {
  for (int dim = 0; dim < 2; dim++) {
    n_[dim] = int(l[dim] / r_cut);
    lc_recip_[dim] = n_[dim] / l[dim];
    if (flag_ext[dim]) {
      origin_[dim] -= l[dim] / n_[dim];
      n_[dim] += 2;
      l_[dim] += 2 * l[dim] / n_[dim];
    }
  }
  n_cells_ = n_.x * n_.y;
  set_padding();
}

CellListBase_2::CellListBase_2(const Vec_2<int> &cell_size,
                               const Vec2d& lc,
                               const Vec2d &gl_l,
                               const Vec2d& origin,
                               const Vec_2<bool>& flag_ext)
                               : n_(cell_size), origin_(origin), lc_recip_(),
                               l_(lc * n_), gl_l_(gl_l), flag_ext_(flag_ext) {
  for (int dim = 0; dim < 2; dim++) {
    lc_recip_[dim] = 1 / lc[dim];
    if (flag_ext[dim]) {
      origin_[dim] -= lc[dim];
      n_[dim] += 2;
      l_[dim] += 2. * lc[dim];
    }
  }
  n_cells_ = n_.x * n_.y;
  set_padding();
}

void CellListBase_2::set_padding() {
  if (flag_ext_.x) {
    real_cell_beg_.x = 1;
    real_cell_end_.x = n_.x - 1;
  } else {
    real_cell_beg_.x = 0;
    real_cell_end_.x = n_.x;
  }
  if (flag_ext_.y) {
    real_cell_beg_.y = 1;
    real_cell_end_.y = n_.y - 1;
  } else {
    real_cell_beg_.y = 0;
    real_cell_end_.y = n_.y;
  }
}

void CellListBase_2::partition(const Vec_2<double>& l, double r_cut,
                               Vec_2<int>& cells_size,
                               Vec_2<double>& cell_len) {
  for (int dim = 0; dim < 2; dim++) {
    cells_size[dim] = int(l[dim] / r_cut);
    cell_len[dim] = l[dim] / cells_size[dim];
  }
}

/**
 * @brief Get a vector to offset the periodic boundary condition
 * !!! The function cannot give correct results for unknown reasons.
 *
 * @param pos             Position of a particle
 * @return Vec_2<double>
 */
Vec_2<double> CellListBase_2::get_offset(const Vec2d& pos) const {
  Vec_2<double> offset{};
  Vec_2<double> dR = pos - origin_;

  //! If canceling the annotation of the following line, the function can give
  //! right results.
  // std::cout << dR << std::endl;

  for (int dim = 0; dim < 2; dim++) {
    if (dR[dim] < 0) {
      offset[dim] = gl_l_[dim];
    } else if (dR[dim] > l_[dim]) {
      offset[dim] = -gl_l_[dim];
    }
  }
  return offset;
}
