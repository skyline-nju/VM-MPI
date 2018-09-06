/**
 * @brief 2d linked-cell list
 *
 * @file cellList2D.h
 * @author skyline-nju
 * @date 2018-9-5
 */
#pragma once
#include <vector>
#include "vect.h"
#include "node.h"

class CellListBase_2 {
public:
  typedef Vec_2<double> Vec2d;

  CellListBase_2(const Vec2d &l,
    double r_cut,
    const Vec2d &gl_l,
    const Vec2d &origin,
    const Vec_2<bool> &flag_ext);

  CellListBase_2(const Vec_2<int> &cell_size,
    const Vec2d &lc,
    const Vec2d &gl_l,
    const Vec2d &origin,
    const Vec_2<bool> &flag_ext);

  static void partition(const Vec_2<double> &l, double r_cut,
    Vec_2<int> &cells_size, Vec_2<double> &cell_len);

  int get_nx(double x) const {
    return int((x - origin_.x) * inverse_lc_.x);
  }
  int get_ny(double y) const {
    return int((y - origin_.y) * inverse_lc_.y);
  }

  template <typename TPar>
  int get_ic(const TPar &p) const {
    return get_nx(p.pos.x) + get_ny(p.pos.y) * n_.x;
  }


  int ncells() const { return ncells_; }
  const Vec2d& origin() const { return origin_; }
  const Vec2d& l() const { return l_; }
  const Vec2d& gl_l() const { return gl_l_; }
  const Vec_2<bool>& flag_ext() const { return flag_ext_; }
  const Vec_2<int>& cells_size() const { return n_; }

  Vec_2<double> get_offset(const Vec2d &pos) const;
protected:
  int ncells_;
  Vec_2<int> n_;
  Vec_2<double> origin_;
  Vec_2<double> inverse_lc_;
  Vec_2<double> l_;
  Vec_2<double> gl_l_;
  Vec_2<bool> flag_ext_;

  const static int cell_offset[4][2];
};

template <typename TNode>
class CellListNode_2 : public CellListBase_2 {
public:
  typedef typename std::vector<TNode*>::iterator IT;
  typedef typename std::vector<TNode*>::const_iterator CIT;
  typedef Vec_2<int> Vec2i;

  CellListNode_2(const Vec2d &l, double r_cut,
                 const Vec2d &gl_l,
                 const Vec2d &origin = Vec2d(),
                 const Vec_2<bool> &flag_ext = Vec_2<bool>())
    : CellListBase_2(l, r_cut, gl_l, origin, flag_ext), head_(ncells_) {
  }

  CellListNode_2(const Vec_2<int> &cell_size,
                 const Vec2d &lc,
                 const Vec2d &gl_l,
                 const Vec2d &origin = Vec2d(),
                 const Vec_2<bool> &flag_ext = Vec_2<bool>())
    : CellListBase_2(cell_size, lc, gl_l, origin, flag_ext), head_(ncells_) {
  }


  template <typename BiFunc1, typename BiFunc2>
  void for_each_pair(BiFunc1 f1, BiFunc2 f2, const Vec2i &ic_beg, const Vec2i &ic_end) const;

  template <typename BiFunc, typename TriFunc>
  void for_each_pair_fast(BiFunc f1, TriFunc f2, const Vec2i &ic_beg, const Vec2i &ic_end) const;

  template <typename BiFunc>
  void for_each_pair_slow(BiFunc f_ij, const Vec2i &ic_beg, const Vec2i &ic_end) const;

  template <typename BiFunc1, typename BiFunc2>
  void for_each_pair(BiFunc1 f1, BiFunc2 f2) const;

  template <typename BiFunc, typename TriFunc>
  void for_each_pair_fast(BiFunc f1, TriFunc f2) const;

  template <typename BiFunc>
  void for_each_pair_slow(BiFunc f_ij) const;

  void create(std::vector<TNode> &p_arr);

  template <typename T>
  void create(std::vector<TNode> &p_arr, T* par_num_arr);

  void recreate(std::vector<TNode> &p_arr);

  template <typename T>
  void recreate(std::vector<TNode> &p_arr, T* par_num_arr);

  template <typename UniFunc>
  void for_each_cell(UniFunc f, const Vec2i &beg, const Vec2i &end);

  void add_node(TNode &p) {
    const int idx = get_ic(p);
    p.append_at_front(&head_[idx]);
  }

  void clear(const Vec2i &first, const Vec2i &last);

  template <typename TPar>
  void replace(BiNode<TPar> &p_new, const BiNode<TPar> &p_old);

  int get_par_num(const Vec2i &beg, const Vec2i &end) const;

  int get_par_num() const;

  void reserve_particles(std::vector<TNode> &p_arr, int new_size, double magnification = 1.1);
protected:
  std::vector<TNode*> head_;
};

template <typename TNode>
template <typename BiFunc>
void CellListNode_2<TNode>::for_each_pair_slow(BiFunc f_ij,
                                              const Vec2i& ic_beg,
                                              const Vec2i& ic_end) const {
  for (int y0 = ic_beg.y; y0 < ic_end.y; y0++) {
    for (int x0 = ic_beg.x; x0 < ic_end.x; x0++) {
      int i0 = x0 + y0 * n_.x;
      if (head_[i0]) {
        for_each_node_pair(head_[i0], f_ij);
        for (int j = 0; j < 4; j++) {
          int y1 = y0 + cell_offset[j][0];
          int x1 = x0 + cell_offset[j][1];
          if (y1 >= n_.y) {
            y1 = 0;
          }
          if (x1 >= n_.x) {
            x1 = 0;
          } else if (x1 < 0) {
            x1 += n_.x;
          }
          int i1 = x1 + y1 * n_.x;
          if (head_[i1]) {
            for_each_node_pair(head_[i0], head_[i1], f_ij);
          }
        }
      }
    }
  }
}

template <typename TNode>
template <typename BiFunc1, typename BiFunc2>
void CellListNode_2<TNode>::for_each_pair(BiFunc1 f1, BiFunc2 f2, const Vec2i& ic_beg, const Vec2i& ic_end) const {
  int y[2];
  int x[2];
  int i[8];

  for (y[0] = ic_beg.y; y[0] < ic_end.y; y[0]++) {
    y[1] = y[0] + 1;
    if (y[1] >= n_.y)
      y[1] -= n_.y;
    const int y_nx[2] = { y[0] * n_.x, y[1] * n_.x };
    for (x[0] = ic_beg.x; x[0] < ic_end.x; x[0]++) {
      x[1] = x[0] + 1;
      if (x[1] >= n_.x)
        x[1] -= n_.x;
      i[0] = x[0] + y_nx[0];
      i[1] = x[1] + y_nx[0];
      i[2] = x[0] + y_nx[1];
      i[3] = x[1] + y_nx[1];

      if (head_[i[0]]) {
        for_each_node_pair(head_[i[0]], f1);
        if (head_[i[1]]) {
          for_each_node_pair(head_[i[0]], head_[i[1]], f2);
        }
        if (head_[i[2]]) {
          for_each_node_pair(head_[i[0]], head_[i[2]], f2);
        }
        if (head_[i[3]]) {
          for_each_node_pair(head_[i[0]], head_[i[3]], f2);
        }
      }
      if (head_[i[1]] && head_[i[2]]) {
        for_each_node_pair(head_[i[1]], head_[i[2]], f2);
      }
    }
  }
}

template <typename TNode>
template <typename BiFunc, typename TriFunc>
void CellListNode_2<TNode>::for_each_pair_fast(BiFunc f1, TriFunc f2,
  const Vec2i& ic_beg, const Vec2i& ic_end) const {
  int y[2];
  int x[2];
  int i[8];

  for (y[0] = ic_beg.y; y[0] < ic_end.y; y[0]++) {
    double ly = 0;
    y[1] = y[0] + 1;
    if (y[1] >= n_.y) {
      y[1] -= n_.y;
      ly = gl_l_.y;
    }
    const int y_nx[2] = { y[0] * n_.x, y[1] * n_.x };
    for (x[0] = ic_beg.x; x[0] < ic_end.x; x[0]++) {
      double lx = 0;
      x[1] = x[0] + 1;
      if (x[1] >= n_.x) {
        x[1] -= n_.x;
        lx = gl_l_.x;
      }
      i[0] = x[0] + y_nx[0];
      i[1] = x[1] + y_nx[0];
      i[2] = x[0] + y_nx[1];
      i[3] = x[1] + y_nx[1];

      if (head_[i[0]]) {
        for_each_node_pair(head_[i[0]], f1);
        if (head_[i[1]]) {
          for_each_node_pair(head_[i[0]], head_[i[1]], Vec2d(lx, 0.), f2);
        }
        if (head_[i[2]]) {
          for_each_node_pair(head_[i[0]], head_[i[2]], Vec2d(0., ly), f2);
        }
        if (head_[i[3]]) {
          for_each_node_pair(head_[i[0]], head_[i[3]], Vec2d(lx, ly), f2);
        }
      }
      if (head_[i[1]] && head_[i[2]]) {
        for_each_node_pair(head_[i[1]], head_[i[2]], Vec2d(-lx, ly), f2);
      }
    }
  }
}

template <typename TNode>
template <typename BiFunc1, typename BiFunc2>
void CellListNode_2<TNode>::for_each_pair(BiFunc1 f1, BiFunc2 f2) const {
  Vec_2<int> beg(0, 0);
  Vec_2<int> end(n_);
  if (flag_ext_.x)
    end.x -= 1;
  if (flag_ext_.y)
    end.y -= 1;
  for_each_pair(f1, f2, beg, end);
}

template <typename TNode>
template <typename BiFunc, typename TriFunc>
void CellListNode_2<TNode>::for_each_pair_fast(BiFunc f1, TriFunc f2) const {
  Vec_2<int> end(n_);
  if (flag_ext_.x)
    end.x -= 1;
  if (flag_ext_.y)
    end.y -= 1;

  for_each_pair_fast(f1, f2, Vec_2<int>(), end);
}

template <typename TNode>
template <typename BiFunc>
void CellListNode_2<TNode>::for_each_pair_slow(BiFunc f_ij) const {
  Vec_2<int> end(n_);
  if (flag_ext_.x)
    end.x -= 1;
  if (flag_ext_.y)
    end.y -= 1;

  for_each_pair_slow(f_ij, Vec_2<int>(), end);
}

template<typename TNode>
void CellListNode_2<TNode>::create(std::vector<TNode>& p_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic = get_ic(*it);
    (*it).append_at_front(&head_[ic]);
  }
}

template <typename TNode>
template <typename T>
void CellListNode_2<TNode>::create(std::vector<TNode>& p_arr, T* par_num_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic = get_ic(*it);
    (*it).append_at_front(&head_[ic]);
    par_num_arr[ic] += 1;
  }
}

template <typename TNode>
void CellListNode_2<TNode>::recreate(std::vector<TNode>& p_arr) {
  for (int ic = 0; ic < ncells_; ic++) {
    head_[ic] = nullptr;
  }
  create(p_arr);
}

template <typename TNode>
template <typename T>
void CellListNode_2<TNode>::recreate(std::vector<TNode>& p_arr, T* par_num_arr) {
  for (int ic = 0; ic < ncells_; ic++) {
    head_[ic] = nullptr;
    par_num_arr[ic] = 0;
  }
  create(p_arr);
}

template <typename TNode>
template <typename UniFunc>
void CellListNode_2<TNode>::for_each_cell(UniFunc f, const Vec2i& beg, const Vec2i& end) {
  for (int y = beg.y; y < end.y; y++) {
    const auto y_nx = y * n_.x;
    for (int x = beg.x; x < end.x; x++) {
      int ic = x + y_nx;
      f(&head_[ic]);
    }
  }
}

template <typename TNode>
void CellListNode_2<TNode>::clear(const Vec2i& first, const Vec2i& last) {
  for (int y0 = first.y; y0 < last.y; y0++) {
    const auto y0_nx = y0 * n_.x;
    for (int x0 = first.x; x0 < last.x; x0++) {
      int ic = x0 + y0_nx;
      head_[ic] = nullptr;
    }
  }
}

template <typename TNode>
template <typename TPar>
void CellListNode_2<TNode>::replace(BiNode<TPar>& p_new, const BiNode<TPar>& p_old) {
  p_new = p_old;
  if (p_new.next) {
    p_new.next->prev = &p_new;
  }
  if (p_new.prev) {
    p_new.prev->next = &p_new;
  } else {
    head_[get_ic(p_new)] = &p_new;
  }
}

template <typename TNode>
int CellListNode_2<TNode>::get_par_num(const Vec2i& beg, const Vec2i& end) const {
  int count = 0;
  for (int y0 = beg.y; y0 < end.y; y0++) {
    const auto y0_nx = y0 * n_.x;
    for (int x0 = beg.x; x0 < end.x; x0++) {
      int ic = x0 + y0_nx;
      TNode * cur_node = head_[ic];
      while (cur_node) {
        count++;
        cur_node = cur_node->next;
      }
    }
  }
  return count;
}

template <typename TNode>
int CellListNode_2<TNode>::get_par_num() const {
  return get_par_num(Vec2i(), n_);
}

template <typename TNode>
void CellListNode_2<TNode>::reserve_particles(std::vector<TNode>& p_arr, int new_size, double magnification) {
  size_t new_cap = new_size * magnification;
  std::cout << "old capacity = " << p_arr.capacity() << "\t";
  p_arr.reserve(new_cap);
  recreate(p_arr);
  std::cout << "new capacity = " << new_cap << std::endl;
}

/**
 * @brief Get a vector that offsets the periodic boundary condition
 *
 * This function works well, however, CellListBase_2::get_offset fails to give
 * expected results, for unknown reasons.
 *
 * @tparam TCellList     Temlate for cell list
 * @param pos            Positon of a particle
 * @param cl             Cell list
 * @return Vec_2<double>
 */
template <typename TCellList>
Vec_2<double> get_offset(const Vec_2<double> &pos, const TCellList &cl) {
  Vec_2<double> offset{};
  Vec_2<double> dR = pos - cl.origin();
  for (int dim = 0; dim < 2; dim++) {
    if (dR[dim] < 0) {
      offset[dim] = cl.gl_l()[dim];
    } else if (dR[dim] > cl.l()[dim]) {
      offset[dim] = -cl.gl_l()[dim];
    }
  }
  return offset;
}