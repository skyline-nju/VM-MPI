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

/**
 * @brief Get a vector that offsets the periodic boundary condition
 *
 * @tparam TCellList     Template for cell list
 * @param pos            Position of a particle
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

  void set_padding();

protected:
  int n_cells_;
  Vec_2<int> n_;
  Vec_2<double> origin_;
  Vec_2<double> lc_recip_;  //  the reciprocal of length of each cell
  Vec_2<double> l_;
  Vec_2<double> gl_l_;
  Vec_2<bool> flag_ext_;

  Vec_2<int> real_cell_beg_{};
  Vec_2<int> real_cell_end_{};

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
    : CellListBase_2(l, r_cut, gl_l, origin, flag_ext), head_(n_cells_) {
  }

  CellListNode_2(const Vec_2<int> &cell_size,
                 const Vec2d &lc,
                 const Vec2d &gl_l,
                 const Vec2d &origin = Vec2d(),
                 const Vec_2<bool> &flag_ext = Vec_2<bool>())
    : CellListBase_2(cell_size, lc, gl_l, origin, flag_ext), head_(n_cells_) {
  }


  template <typename BiFunc1, typename BiFunc2>
  void for_each_pair(BiFunc1 f1, BiFunc2 f2, const Vec2i &ic_beg, const Vec2i &ic_end) const;
  template <typename BiFunc1, typename BiFunc2>
  void for_each_pair(BiFunc1 f1, BiFunc2 f2) const { for_each_pair(f1, f2, Vec_2<int>(), real_cell_end_); }

  template <typename BiFunc, typename TriFunc>
  void for_each_pair_fast(BiFunc f1, TriFunc f2, const Vec2i &ic_beg, const Vec2i &ic_end) const;
  template <typename BiFunc, typename TriFunc>
  void for_each_pair_fast(BiFunc f1, TriFunc f2) const { for_each_pair_fast(f1, f2, Vec_2<int>(), real_cell_end_); }

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
  int get_par_num() const;
  template <typename T>
  int get_par_num(const std::vector<T> &n_arr) const;

  void reserve_particles(std::vector<TNode> &p_arr, int new_size, double magnification = 1.1);

  template <typename UniFunc, typename TRan, typename T1, typename T2>
  void birth_death(std::vector<TNode> &p_arr, UniFunc f_rate, TRan &myran,
                   std::vector<T1> &n_arr, std::vector<Vec_2<T2>> &v_arr);

  template <typename UniFunc, typename TRan, typename T1, typename T2>
  void birth_death(std::vector<TNode> &p_arr, UniFunc f_rate, TRan &myran,
                   std::vector<T1> &n_arr, std::vector<Vec_2<T2>> &v_arr,
                   const Vec_2<int> &cg_box_l);

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
  Vec_2<int> beg = real_cell_beg_;
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
  //! sorted in descending order
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
int CellListNode_2<TNode>::get_par_num() const {
  return get_par_num(Vec2i(), n_);
}

template <typename TNode>
template <typename T>
int CellListNode_2<TNode>::get_par_num(const std::vector<T>& n_arr) const {
  int n_sum = 0;
  for (int row = real_cell_beg_.y; row < real_cell_end_.y; row++) {
    for (int col = real_cell_beg_.x; col < real_cell_end_.x; col++) {
      int ic = col + row * n_.x;
      n_sum += n_arr[ic];
    }
  }
  return n_sum;
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



template<typename TNode>
template <typename UniFunc, typename TRan, typename T1, typename T2>
void CellListNode_2<TNode>::birth_death(std::vector<TNode> &p_arr, UniFunc f_rate, TRan &myran,
                                        std::vector<T1> &n_arr, std::vector<Vec_2<T2>> &v_arr) {
  std::vector<int> vacant_pos;
  vacant_pos.reserve(2048);
  for (int row = real_cell_beg_.y; row < real_cell_end_.y; row++) {
    const int row_nx = row * n_.x;
    for (int col = real_cell_beg_.x; col < real_cell_end_.x; col++) {
      const int ic = col + row_nx;
      const int my_n = n_arr[ic];
      double rate = f_rate(double(my_n));
      if (rate < 0.) {
        if (myran.doub() < -rate) {
          const int k = int(myran.doub() * my_n);
          vacant_pos.push_back(remove_particle(p_arr, ic, k));
        }
      } else if (rate > 0.) {
        if (myran.doub() < rate) {
          Vec_2<double> new_pos(col + myran.doub(), row + myran.doub());
          new_pos += origin_;
          Vec_2<double> new_v{};
          if (my_n == 0) {
            const double theta = PI * 2. * myran.doub();
            new_v.x = std::cos(theta);
            new_v.y = std::sin(theta);
          } else {
            new_v.x = v_arr[ic].x / my_n;
            new_v.y = v_arr[ic].y / my_n;
            new_v.normalize();
          }
          add_particle(p_arr, vacant_pos, new_pos, new_v);
        }
      }
    }
  }
  make_compact(p_arr, vacant_pos);
}

template <typename TNode>
template <typename UniFunc, typename TRan, typename T1, typename T2>
void CellListNode_2<TNode>::birth_death(std::vector<TNode>& p_arr, UniFunc f_rate, TRan& myran,
                                        std::vector<T1>& n_arr, std::vector<Vec_2<T2>>& v_arr,
                                        const Vec_2<int>& cg_box_l) {
  const int cg_nrows = (real_cell_end_.y - real_cell_beg_.y) / cg_box_l.y;
  const int cg_ncols = (real_cell_end_.x - real_cell_beg_.x) / cg_box_l.x;
  const double area = cg_box_l.x * cg_box_l.y;

  std::vector<int> vacant_pos;
  vacant_pos.reserve(2048);
  std::vector<int> n_box(cg_box_l.x * cg_box_l.y);

  for (int cg_row = 0; cg_row < cg_nrows; cg_row++) {
    const int fine_row_beg = cg_row * cg_box_l.y + real_cell_beg_.y;
    const int fine_row_end = fine_row_beg + cg_box_l.y;
    for (int cg_col = 0; cg_col < cg_ncols; cg_col++) {
      const int fine_col_beg = cg_col * cg_box_l.x + real_cell_beg_.x;
      const int fine_col_end = fine_col_beg + cg_box_l.x;
      int cg_num = 0;
      Vec_2<double> cg_v{};
      int idx_box = 0;
      for (int fine_row = fine_row_beg; fine_row < fine_row_end; fine_row++) {
        const int fine_row_nx = fine_row * n_.x;
        for (int fine_col = fine_col_beg; fine_col <fine_col_end; fine_col++) {
          const int ic = fine_col + fine_row_nx;
          cg_num += n_arr[ic];
          if (cg_v.dot(v_arr[ic]) >= 0.) {
            cg_v += v_arr[ic];
          } else {
            cg_v  -= v_arr[ic];
          }
          n_box[idx_box] = n_arr[ic];
          idx_box++;
        }
      }
      if (cg_num > 0) {
        double cg_rho = cg_num / area;
        double rate = f_rate(cg_rho) * cg_num;
        if (rate < 0.) {
          do {
            if (rate <= -1. || myran.doub() < -rate) {
              int k = cg_num * myran.doub();
              int i = 0;
              while (k >= n_box[i]) {
                k -= n_box[i];
                i++;
              }
              int ic = i % cg_box_l.x + fine_col_beg + n_.x * (i / cg_box_l.x + fine_row_beg);
              vacant_pos.push_back(remove_particle(p_arr, ic, k));
              n_box[i]--;
            }
            rate += 1.;
          } while (rate < 0.);
          
        } else {
          do {
            if (rate >= 1. || myran.doub() < rate) {
              Vec_2<double> new_pos(fine_col_beg + myran.doub() * cg_box_l.x,
                fine_row_beg + myran.doub() * cg_box_l.y);
              new_pos += origin_;
              Vec_2<double> new_v{ cg_v.x / cg_num, cg_v.y / cg_num };
              new_v.normalize();
              add_particle(p_arr, vacant_pos, new_pos, new_v);
            }
            rate -= 1;
          } while (rate > 0.);
          
        }
      }
    }
  }
  make_compact(p_arr, vacant_pos);
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
