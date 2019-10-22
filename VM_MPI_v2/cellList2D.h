/**
 * @brief 2d linked-cell list
 *
 * @file cellList2D.h
 * @author skyline-nju
 * @date 2018-9-5
 */
#pragma once
#include <vector>
#include <algorithm>
#include "vect.h"
#include "node.h"
#include "comn.h"

template <typename T>
struct RectBlock_2 {
  RectBlock_2() = default;

  RectBlock_2(const Vec_2<T>& n, const Vec_2<bool>& need_pad, const T pad_w);

  T size_x() const { return end.x - beg.x; }
  T size_y() const { return end.y - beg.y; }

  Vec_2<T> beg;
  Vec_2<T> end;
};

template <typename T>
RectBlock_2<T>::RectBlock_2(const Vec_2<T>& n, const Vec_2<bool>& need_pad, const T pad_w) {
  if (need_pad.x) {
    beg.x = pad_w;
    end.x = n.x - pad_w;
  } else {
    beg.x = 0;
    end.x = n.x;
  }
  if (need_pad.y) {
    beg.y = pad_w;
    end.y = n.y - pad_w;
  } else {
    beg.y = 0;
    end.y = n.y;
  }
}

class CellListBase_2 {
public:
  typedef Vec_2<double> Vec2d;

  CellListBase_2(const Vec_2<int> &cell_size,
                 const Vec2d &lc,
                 const Vec2d &gl_l,
                 const Vec2d &origin,
                 const Vec_2<bool> &flag_ext);

  int get_nx(double x) const { return int((x - origin_.x) * lc_recip_.x); }
  int get_ny(double y) const { return int((y - origin_.y) * lc_recip_.y); }
  template <typename TPar>
  int get_ic(const TPar &p) const { return get_nx(p.pos.x) + get_ny(p.pos.y) * n_.x; }

  int n_cells() const { return n_cells_; }
  const Vec2d& origin() const { return origin_; }
  const Vec2d& l() const { return l_; }
  const Vec2d& gl_l() const { return gl_l_; }
  const Vec_2<bool>& flag_ext() const { return flag_ext_; }
  const Vec_2<int>& cells_size() const { return n_; }

  template <typename T>
  int get_par_num(const std::vector<T>& n_arr) const;

protected:
  int n_cells_;
  Vec_2<int> n_;
  Vec_2<double> origin_;
  Vec_2<double> lc_recip_;  //  the reciprocal of length of one cell
  Vec_2<double> l_;
  Vec_2<double> gl_l_;
  Vec_2<bool> flag_ext_;

  RectBlock_2<int> real_block_;

  const int cell_offset[4][2] = { { 0,  1 },
                                  { 1, -1 },
                                  { 1,  0 },
                                  { 1,  1 }};
};


inline CellListBase_2::CellListBase_2(const Vec_2<int>& cell_size,
                                      const Vec2d& lc,
                                      const Vec2d& gl_l,
                                      const Vec2d& origin,
                                      const Vec_2<bool>& flag_ext)
                                      : n_(cell_size), origin_(origin), lc_recip_(1./lc.x, 1./lc.y),
                                      l_(lc* n_), gl_l_(gl_l), flag_ext_(flag_ext),
                                      real_block_(cell_size, flag_ext_, 1) {
  for (int dim = 0; dim < 2; dim++) {
    if (flag_ext[dim]) {
      origin_[dim] -= lc[dim];
      n_[dim] += 2;
      l_[dim] += 2. * lc[dim];
    }
  }
  n_cells_ = n_.x * n_.y;
}

template <typename T>
int CellListBase_2::get_par_num(const std::vector<T>& n_arr) const {
  int n_sum = 0;
  for (int row = real_block_.beg.y; row < real_block_.end.y; row++) {
    for (int col = real_block_.beg.x; col < real_block_.end.x; col++) {
      int ic = col + row * n_.x;
      n_sum += n_arr[ic];
    }
  }
  return n_sum;
}

template <typename TNode>
class CellListNode_2 : public CellListBase_2 {
public:
  typedef typename std::vector<TNode*>::iterator IT;
  typedef typename std::vector<TNode*>::const_iterator CIT;
  typedef Vec_2<int> Vec2i;

  CellListNode_2(const Vec_2<int> &cell_size,
                 const Vec2d &lc,
                 const Vec2d &gl_l,
                 const Vec2d &origin = Vec2d(),
                 const Vec_2<bool> &flag_ext = Vec_2<bool>())
    : CellListBase_2(cell_size, lc, gl_l, origin, flag_ext), head_(n_cells_) {
  }
  
  template <typename TDomain>
  CellListNode_2(const TDomain& dm) : CellListBase_2(dm.grid().n(), dm.grid().lc(), dm.gl_l(), dm.origin(), dm.flag_comm()),
                                      head_(n_cells_) {}

  template <typename BiFunc1, typename BiFunc2>
  void for_each_pair(BiFunc1 f1, BiFunc2 f2, const Vec2i &ic_beg, const Vec2i &ic_end) const;
  template <typename BiFunc1, typename BiFunc2>
  void for_each_pair(BiFunc1 f1, BiFunc2 f2) const { for_each_pair(f1, f2, Vec_2<int>(), real_block_.end); }

  template <typename BiFunc, typename TriFunc>
  void for_each_pair_fast(BiFunc f1, TriFunc f2, const Vec2i &ic_beg, const Vec2i &ic_end) const;
  template <typename BiFunc, typename TriFunc>
  void for_each_pair_fast(BiFunc f1, TriFunc f2) const { for_each_pair_fast(f1, f2, Vec_2<int>(), real_block_.end); }

  template <typename BiFunc1, typename BiFunc2>
  void for_each_pair_slow(BiFunc1 f1, BiFunc2 f2, const Vec2i &ic_beg, const Vec2i &ic_end) const;
  template <typename BiFunc1, typename BiFunc2>
  void for_each_pair_slow(BiFunc1 f1, BiFunc2 f2) const;

  template <typename UniFunc>
  void for_each_cell(UniFunc f, const Vec2i &beg, const Vec2i &end);

  void create(std::vector<TNode> &p_arr);
  void recreate(std::vector<TNode> &p_arr);

  template <typename T1, typename T2>
  void create(std::vector<TNode> &p_arr, std::vector<T1> &num_arr, std::vector<Vec_2<T2>> &v_arr);
  template <typename T1, typename T2>
  void recreate(std::vector<TNode> &p_arr, std::vector<T1> &num_arr, std::vector<Vec_2<T2>> &v_arr);

  void add_node(TNode &p) { p.append_at_front(&head_[get_ic(p)]); }
  template <typename T1, typename T2>
  void add_node(TNode &p, std::vector<T1> &n_arr, std::vector<T2> &v_arr);

  void clear(const Vec2i &first, const Vec2i &last);

  void replace(TNode &p1, const TNode &p2);

  void make_compact(std::vector<TNode>& p_arr, std::vector<int>& vacancy);

  int get_par_num(const Vec2i &beg, const Vec2i &end) const;
  int get_par_num() const { return get_par_num(Vec2i(), n_); }

  void reserve_particles(std::vector<TNode> &p_arr, int new_size, double magnification = 1.1);

  int remove_particle(std::vector<TNode> &p_arr, int ic, int ip);

  void add_particle(std::vector<TNode> &p_arr, std::vector<int> &vacant_pos,
                    const Vec_2<double> &pos, const Vec_2<double> &vel);
protected:
  std::vector<TNode*> head_;
};

template <typename TNode>
template <typename BiFunc1, typename BiFunc2>
void CellListNode_2<TNode>::for_each_pair_slow(BiFunc1 f1, BiFunc2 f2,
                                              const Vec2i& ic_beg,
                                              const Vec2i& ic_end) const {
  for (int y0 = ic_beg.y; y0 < ic_end.y; y0++) {
    for (int x0 = ic_beg.x; x0 < ic_end.x; x0++) {
      int i0 = x0 + y0 * n_.x;
      if (head_[i0]) {
        for_each_node_pair(head_[i0], f1);
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
            for_each_node_pair(head_[i0], head_[i1], f2);
          }
        }
      }
    }
  }
}

template <typename TNode>
template <typename BiFunc1, typename BiFunc2>
void CellListNode_2<TNode>::for_each_pair_slow(BiFunc1 f1, BiFunc2 f2) const {
  Vec_2<int> beg = real_block_.beg;
  Vec_2<int> end = n_;
  if (flag_ext_.x) {
    beg.x = 0;
    end.x = n_.x;
  }
  if (flag_ext_.y) {
    beg.y -= 1;
    end.y -= 1;
  }
  for_each_pair_slow(f1, f2, beg, end);
}

template <typename TNode>
template <typename BiFunc1, typename BiFunc2>
void CellListNode_2<TNode>::for_each_pair(BiFunc1 f1, BiFunc2 f2,
                                          const Vec2i& ic_beg,
                                          const Vec2i& ic_end) const {
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
                                               const Vec2i& ic_beg,
                                               const Vec2i& ic_end) const {
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

template<typename TNode>
void CellListNode_2<TNode>::create(std::vector<TNode>& p_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    add_node(*it);
  }
}

template <typename TNode>
void CellListNode_2<TNode>::recreate(std::vector<TNode>& p_arr) {
  for (int ic = 0; ic < n_cells_; ic++) {
    head_[ic] = nullptr;
  }
  create(p_arr);
}

template <typename TNode>
template <typename T1, typename T2>
void CellListNode_2<TNode>::create(std::vector<TNode>& p_arr,
                                   std::vector<T1>& num_arr,
                                   std::vector<Vec_2<T2>>& v_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    add_node(*it, num_arr, v_arr);
  }
}

template <typename TNode>
template <typename T1, typename T2>
void CellListNode_2<TNode>::recreate(std::vector<TNode>& p_arr,
                                     std::vector<T1>& num_arr,
                                     std::vector<Vec_2<T2>>& v_arr) {
  for (int ic = 0; ic < n_cells_; ic++) {
    head_[ic] = nullptr;
    num_arr[ic] = 0;
    v_arr[ic].x = v_arr[ic].y = 0;
  }
  create(p_arr, num_arr, v_arr);
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
template <typename T1, typename T2>
void CellListNode_2<TNode>::add_node(TNode& p, std::vector<T1>& n_arr, std::vector<T2>& v_arr) {
  const int ic = get_ic(p);
  n_arr[ic] += 1;
#ifdef POLAR_ALIGN
  v_arr[ic] += p.ori;
#else
  if (v_arr[ic].dot(p.ori) >= 0.) {
    v_arr[ic] += p.ori;
  } else {
    v_arr[ic] -= p.ori;
  }
#endif
  p.append_at_front(&head_[ic]);
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
void CellListNode_2<TNode>::replace(TNode& p1, const TNode& p2) {
  p1 = p2;
  if (p1.next) {
    p1.next->prev = &p1;
  }
  if (p1.prev) {
    p1.prev->next = &p1;
  } else {
    head_[get_ic(p1)] = &p1;
  }
}


template<typename TNode>
void CellListNode_2<TNode>::make_compact(std::vector<TNode>& p_arr,
                                         std::vector<int>& vacancy) {
  //! vacancy should be sorted in descending order
  std::sort(vacancy.begin(), vacancy.end(), std::greater<int>());
  int k = 0;
  while (k < vacancy.size()) {
    if (p_arr.size() - 1 == vacancy[k]) {
      p_arr.pop_back();
      k++;
    } else {
      auto i = vacancy.back();
      replace(p_arr[i], p_arr.back());
      p_arr.pop_back();
      vacancy.pop_back();
    }
  }
  vacancy.clear();
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
void CellListNode_2<TNode>::reserve_particles(std::vector<TNode>& p_arr,
                                              int new_size,
                                              double magnification) {
  size_t new_cap = new_size * magnification;
  std::cout << "old capacity = " << p_arr.capacity() << "\t";
  p_arr.reserve(new_cap);
  recreate(p_arr);
  std::cout << "new capacity = " << new_cap << std::endl;
}

template <typename TNode>
int CellListNode_2<TNode>::remove_particle(std::vector<TNode>& p_arr, int ic, int ip) {
  int real_ip = 0;
  TNode *cur_node = head_[ic];
  int count = 0;
  do {
    if (count == ip) {
      if (count == 0) {
        head_[ic] = cur_node->next;
        if (cur_node->next) {
          cur_node->next->prev = nullptr;
        }
      } else {
        cur_node->prev->next = cur_node->next;
        if (cur_node->next) {
          cur_node->next->prev = cur_node->prev;
        }
      }
      real_ip = cur_node - &p_arr[0];
      break;
    }
    cur_node = cur_node->next;
    count++;
  } while (cur_node);
  return real_ip;
}

template <typename TNode>
void CellListNode_2<TNode>::add_particle(std::vector<TNode>& p_arr, std::vector<int>& vacant_pos,
                                         const Vec_2<double>& pos, const Vec_2<double>& vel) {
  int idx;
  if (vacant_pos.empty()) {
    idx = p_arr.size();
    p_arr.push_back(TNode()); // need improve
  } else {
    idx = vacant_pos.back();
    vacant_pos.pop_back();
  }
  p_arr[idx].update(pos, vel);
  add_node(p_arr[idx]);
}

/**
 * @brief Get a vector that offsets the periodic boundary condition
 *
 * @tparam TCellList     Template for cell list
 * @param pos            Position of a particle
 * @param cl             Cell list
 * @return Vec_2<double>
 */
template <typename TCellList>
Vec_2<double> get_offset(const Vec_2<double>& pos, const TCellList& cl) {
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