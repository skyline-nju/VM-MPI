#ifndef CELLLIST2D_H
#define CELLLIST2D_H
#include <list>
#include <vector>
#include <iostream>
#include <typeinfo>
#include "vect.h"
#include "comn.h"
#include "node.h"
class _CellListBase_2 {
public:
  _CellListBase_2(double Lx, double Ly, double x0 = 0, double y0 = 0,
    double r_cut = 1, bool comm_x = false, bool comm_y = false);

  template <class _TPar>
  int get_ic(const _TPar &p) const;

  template <class _TPar, class BiFunc>
  void for_each_pair0(const std::vector<_TPar> &cl, BiFunc f_ij,
    int row_beg, int row_end, int col_beg, int col_end);

  template <class _TPar>
  void create_lists(std::vector<std::list<_TPar *>> &cl, std::vector<_TPar> &p_arr);

  template <class _TNode>
  void create_lists(std::vector<_TNode* > &cl, std::vector<_TNode> &p_arr);

  template <class _TNode, class UniFunc>
  void update_lists(std::vector<_TNode* > &cl, UniFunc move);

protected:
  Vec_2<double> origin;
  Vec_2<double> l_box;
  Vec_2<double> inverse_lc;
  Vec_2<double> l_cell;
  Vec_2<int> n_bins;
  int ncells;
};

template<class _TPar>
inline int _CellListBase_2::get_ic(const _TPar & p) const {
  return int((p.x - origin.x) * inverse_lc.x)
    + int((p.y - origin.y) * inverse_lc.y) * n_bins.x;
}

template<class _TPar, class BiFunc>
void _CellListBase_2::for_each_pair0(const std::vector<_TPar>& cl, BiFunc f_ij,
                                     int row_beg, int row_end, int col_beg, int col_end) {
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
      if (not_empty(cl[ic0])) {
        for_each_node_pair(cl[ic0], f_ij);
        int ic3 = col_right + row_upper_times_ncols;
        bool flag_c1_c2 = true;
        if (not_empty(cl[ic1]))
          for_each_node_pair(cl[ic0], cl[ic1], f_ij);
        else
          flag_c1_c2 = false;
        if (not_empty(cl[ic2]))
          for_each_node_pair(cl[ic0], cl[ic2], f_ij);
        else
          flag_c1_c2 = false;
        if (flag_c1_c2)
          for_each_node_pair(cl[ic1], cl[ic2], f_ij);
        if (not_empty(cl[ic3]))
          for_each_node_pair(cl[ic0], cl[ic3], f_ij);
      } else if (not_empty(cl[ic1]) && not_empty(cl[ic2])) {
        for_each_node_pair(cl[ic1], cl[ic2], f_ij);
      }
    }
  }
}

template<class _TPar>
void _CellListBase_2::create_lists(std::vector<std::list<_TPar*>>& cl,
                                   std::vector<_TPar>& p_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic = get_ic(*it);
    cl[ic].push_back(&(*it));
  }
}

template <class _TNode>
void _CellListBase_2::create_lists(std::vector<_TNode* > &cl,
                                   std::vector<_TNode> &p_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic = get_ic(*it);
    (*it).append_at_front(&cl[ic]);
  }
}

template <class _TNode, class UniFunc>
void _CellListBase_2::update_lists(std::vector<_TNode* > &cl, UniFunc move) {
  std::vector<_TNode*> tail(ncells, nullptr);
  std::vector<_TNode*> head_tmp(ncells, nullptr);
  for (int ic = 0; ic < ncells; ic++) {
    if (cl[ic]) {
      auto cur_node = cl[ic];
      _TNode * pre_node = nullptr;
      do {
        move(cur_node);
        int ic_new = get_ic(*cur_node);
        if (ic_new != ic) {
          cur_node->break_away(&cl[ic], pre_node);
          auto tmp_node = cur_node;
          cur_node = cur_node->next;
          tmp_node->append_at_front(&head_tmp[ic_new]);
        } else {
          pre_node = cur_node;
          cur_node = cur_node->next;
        }
      } while (cur_node);
      tail[ic] = pre_node;
    }
  }
  for (int ic = 0; ic < ncells; ic++) {
    merge_list(&cl[ic], &tail[ic], &head_tmp[ic]);
  }

}
// use std::list to link particles
template <class _TPar>
class CellList_2 : public _CellListBase_2 {
public:
  CellList_2(double Lx, double Ly, double x0 = 0, double y0 = 0,
    double r_cut = 1, bool comm_x = false, bool comm_y = false);

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij) {
    for_each_pair0(cell, f_ij, 0, n_bins.y, 0, n_bins.x); }

  void create(std::vector<_TPar> &p_arr);

  void recreate(std::vector<_TPar> &p_arr) {
    del_list(cell);
    create(p_arr);
  }

  void update();

  template <class UniFunc>
  void update(UniFunc move);

private:
  std::vector<std::list<_TPar*> > cell;
  std::vector<int> list_len;
};

template<class _TPar>
CellList_2<_TPar>::CellList_2(double Lx, double Ly, double x0, double y0,
                              double r_cut, bool comm_x, bool comm_y):
                              _CellListBase_2(Lx, Ly, x0, y0, r_cut, comm_x, comm_y) {
  cell.reserve(ncells);
  list_len.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.emplace_back();
    list_len.push_back(0);
  }
  std::cout << n_bins.x << "\t" << n_bins.y << "\n";
}

template<class _TPar>
void CellList_2<_TPar>::create(std::vector<_TPar>& p_arr) {
  create_lists(cell, p_arr);
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template<class _TPar>
void CellList_2<_TPar>::update() {
  for (int iy = 0; iy < n_bins.y; iy++) {
    int iy_times_nx = iy * n_bins.x;

    for (int ix = 0; ix < n_bins.x; ix++) {
      int ic = ix + iy_times_nx;
      int depth = list_len[ic];
      if (depth) {
        int it_count = 0;
        for (auto it = cell[ic].begin(); it_count < depth; it_count++) {
          _TPar *ptr = *it;
          int ic_new = get_ic(*ptr);
          if (ic_new == ic) {
            ++it;
          } else {
            it = cell[ic].erase(it);
            cell[ic_new].push_back(ptr);
          }
        }
      }
    }
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template<class _TPar>
template<class UniFunc>
void CellList_2<_TPar>::update(UniFunc move) {
  for (int ic = 0; ic < ncells; ic++) {
    int depth = list_len[ic];
    if (depth) {
      int it_count = 0;
      for (auto it = cell[ic].begin(); it_count < depth; it_count++) {
        _TPar *ptr = *it;
        move(ptr);
        int ic_new = get_ic(*ptr);
        if (ic_new == ic) {
          ++it;
        } else {
          it = cell[ic].erase(it);
          cell[ic_new].push_back(ptr);
        }
      }
    }
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

// wrap the particle with class UniNode
template <class _TPar>
class CellList_2<UniNode<_TPar>> : public _CellListBase_2 {
public:
  CellList_2(double Lx, double Ly, double x0 = 0, double y0 = 0,
    double r_cut = 1, bool comm_x = false, bool comm_y = false);

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij) {
    for_each_pair0(cell, f_ij, 0, n_bins.y, 0, n_bins.x); }

  void create(std::vector<UniNode<_TPar>> &p_arr) { create_lists(cell, p_arr); }

  void recreate(std::vector<UniNode<_TPar>> &p_arr) {
    del_list(cell);
    create(p_arr);
  }

  template <class UniFunc>
  void update(UniFunc move) {update_lists(cell, move); }

private:
  std::vector<UniNode<_TPar> *> cell;
};

template <class _TPar>
CellList_2<UniNode<_TPar>>::CellList_2(double Lx, double Ly, double x0, double y0,
                                       double r_cut, bool comm_x, bool comm_y):
  _CellListBase_2(Lx, Ly, x0, y0, r_cut, comm_x, comm_y) {
  cell.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.emplace_back();
  }
  std::cout << n_bins.x << "\t" << n_bins.y << "\n";
}

// wrap the particle with BNode
template <class _TPar>
class CellList_2<BNode<_TPar>> : public _CellListBase_2 {
public:
  CellList_2(double Lx, double Ly, double x0 = 0, double y0 = 0,
    double r_cut = 1, bool comm_x = false, bool comm_y = false);

  bool has_member(int ic) const { return cell[ic] != nullptr; }

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij) {
    for_each_pair0(cell, f_ij, 0, n_bins.y, 0, n_bins.x); }

  void create(std::vector<BNode<_TPar>> &p_arr) { create_lists(cell, p_arr); }

  void recreate(std::vector<BNode<_TPar>> &p_arr) {
    del_list(cell);
    create(p_arr);
  }

  template <class UniFunc>
  void update(UniFunc move) { update_lists(cell, move); }

private:
  std::vector<BNode<_TPar> *> cell;
};

template <class _TPar>
CellList_2<BNode<_TPar>>::CellList_2(double Lx, double Ly, double x0, double y0,
                                     double r_cut, bool comm_x, bool comm_y):
  _CellListBase_2(Lx, Ly, x0, y0, r_cut, comm_x, comm_y) {
  cell.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.emplace_back();
  }
  std::cout << n_bins.x << "\t" << n_bins.y << "\n";
}

#endif