#ifndef CELLLIST2D_H
#define CELLLIST2D_H
#include <list>
#include <vector>
#include <iterator>
#include <iostream>
#include <typeinfo>
#include "vect.h"
#include "comn.h"
#include "node.h"

template <typename _TNode>
class Cell_list_2 {
public:
  typedef typename std::vector<_TNode*>::iterator IT;
  typedef typename std::vector<_TNode*>::const_iterator CIT;
  Cell_list_2(double Lx, double Ly, double x0 = 0, double y0 = 0, double r_cut = 1,
    bool comm_x = false, bool comm_y = false);

  int get_ic(const _TNode &p) const;

  template <typename BiFunc>
  void for_each_pair(BiFunc f_ij, int row_beg = 0, int row_end = -1, int col_beg = 0, int col_end = -1);

  void create(std::vector<_TNode> &p_arr);

  void recreate(std::vector<_TNode> &p_arr);


  template <typename UniFunc>
  void update(UniFunc move, int first_ic, CIT first, CIT last);

  template <typename UniFunc>
  void update(UniFunc move, int ic, _TNode* head0);

  template <typename UniFunc>
  void update(UniFunc move, int first_ic, const std::vector<_TNode*> &head0) {
    update(move, first_ic, head0.cbegin(), head0.cend());
  }

  template <typename UniFunc>
  void update(UniFunc move);

  template <typename UniFunc>
  void update_by_row(UniFunc move);

protected:
  Vec_2<double> origin;
  Vec_2<double> l_box;
  Vec_2<double> inverse_lc;
  Vec_2<double> l_cell;
  Vec_2<int> n_bins;
  int ncells;
  std::vector<_TNode*> head;
};

template<typename _TNode>
Cell_list_2<_TNode>::Cell_list_2(double Lx, double Ly, double x0, double y0, double r_cut,
  bool comm_x, bool comm_y) {
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
  head.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    head.emplace_back();
  }
}

template<typename _TNode>
inline int Cell_list_2<_TNode>::get_ic(const _TNode & p) const {
  return int((p.x - origin.x) * inverse_lc.x)
    + int((p.y - origin.y) * inverse_lc.y) * n_bins.x;
}

template<typename _TNode>
void Cell_list_2<_TNode>::create(std::vector<_TNode>& p_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic = get_ic(*it);
    (*it).append_at_front(&head[ic]);
  }
}

template<typename _TNode>
void Cell_list_2<_TNode>::recreate(std::vector<_TNode>& p_arr) {
  for (int ic = 0; ic < ncells; ic++) {
    head[ic] = nullptr;
  }
  create(p_arr);
}

template<typename _TNode>
template<typename BiFunc>
void Cell_list_2<_TNode>::for_each_pair(BiFunc f_ij, int row_beg, int row_end,
  int col_beg, int col_end) {
  if (row_end == -1)
    row_end = n_bins.y;
  if (col_end == -1)
    col_end = n_bins.x;
  for (int row = row_beg; row < row_end; row++) {
    int row_upper = row + 1;
    if (row_upper == n_bins.y)
      row_upper = 0;
    int row_times_ncols = row * n_bins.x;
    int row_upper_times_ncols = row_upper * n_bins.x;
    for (int col = col_beg; col < col_end; col++) {
      int ic0 = col + row_times_ncols;
      int col_right = col + 1;
      if (col_right == n_bins.x)
        col_right = 0;
      int ic1 = col_right + row_times_ncols;
      int ic2 = col + row_upper_times_ncols;
      if (head[ic0]) {
        for_each_node_pair(head[ic0], f_ij);
        int ic3 = col_right + row_upper_times_ncols;
        bool flag_c1_c2 = true;
        if (head[ic1])
          for_each_node_pair(head[ic0], head[ic1], f_ij);
        else
          flag_c1_c2 = false;
        if (head[ic2])
          for_each_node_pair(head[ic0], head[ic2], f_ij);
        else
          flag_c1_c2 = false;
        if (flag_c1_c2)
          for_each_node_pair(head[ic1], head[ic2], f_ij);
        if (head[ic3])
          for_each_node_pair(head[ic0], head[ic3], f_ij);
      } else if (head[ic1] && head[ic2]) {
        for_each_node_pair(head[ic1], head[ic2], f_ij);
      }
    }
  }
}

template<typename _TNode>
template<typename UniFunc>
void Cell_list_2<_TNode>::update_by_row(UniFunc move) {
  int nx = n_bins.x;
  std::vector<_TNode*> pre_row(head.cend() - nx, head.cend());
  auto first = head.cbegin();
  auto last = first + nx;
  std::vector<_TNode*> cur_row(first, last);
  first += nx;
  last += nx;
  std::vector<_TNode*> nxt_row(first, last);
  update(move, 0, cur_row);
  cur_row.swap(nxt_row);

  for (int row = 1; row < n_bins.y - 2; row++) {
    first += nx;
    last += nx;
    nxt_row.assign(first, last);
    update(move, row * nx, cur_row.cbegin(), cur_row.cend());
    cur_row.swap(nxt_row);
  }
  update(move, (n_bins.y - 2) * nx, cur_row);
  update(move, (n_bins.y - 1) * nx, pre_row);
}

template <typename _TNode>
template <typename UniFunc>
void Cell_list_2<_TNode>::update(UniFunc move, int ic, _TNode* head0) {
  auto cur_node = head0;
  bool flag = cur_node == head[ic] ? true : false;
  do {
    move(cur_node);
    int ic_new = get_ic(*cur_node);
    if (ic_new != ic) {
      if (flag) {
        cur_node->break_away(&head[ic]);
      } else {
        cur_node->break_away();
      }
      auto tmp_node = cur_node;
      cur_node = cur_node->next;
      tmp_node->append_at_front(&head[ic_new]);
    } else {
      cur_node = cur_node->next;
    }
  } while (cur_node);
}

template<typename _TNode>
template<typename UniFunc>
void Cell_list_2<_TNode>::update(UniFunc move, int first_ic, CIT first, CIT last) {
  int ic = first_ic;
  for (auto it = first; it != last; ++it) {
    if (*it) {
      update(move, ic, *it);
    }
    ic++;
  }
}

template<typename _TNode>
template<typename UniFunc>
inline void Cell_list_2<_TNode>::update(UniFunc move) {
  std::vector<_TNode*> head_cur(head);
  update(move, 0, head_cur.cbegin(), head_cur.cend());
}
#endif

