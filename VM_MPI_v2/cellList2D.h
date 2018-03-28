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

  virtual bool has_member(int ic) const = 0;

  template <class UniFunc, class BiFunc>
  void for_each_nearest_cell(UniFunc intra_cell, BiFunc inter_cell,
    int row_beg, int row_end, int col_beg, int col_end);

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

template<class UniFunc, class BiFunc>
void _CellListBase_2::for_each_nearest_cell(UniFunc intra_cell, BiFunc inter_cell,
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
      if (has_member(ic0)) {
        intra_cell(ic0);
        int ic3 = col_right + row_upper_times_ncols;
        bool flag_c1_c2 = true;
        if (has_member(ic1))
          inter_cell(ic0, ic1);
        else
          flag_c1_c2 = false;
        if (has_member(ic2))
          inter_cell(ic0, ic2);
        else
          flag_c1_c2 = false;
        if (flag_c1_c2)
          inter_cell(ic1, ic2);
        if (has_member(ic3))
          inter_cell(ic0, ic3);
      } else if (has_member(ic1) && has_member(ic2)) {
        inter_cell(ic1, ic2);
      }   
    }
  }
}

template <class _TPar>
class CellList_2 : public _CellListBase_2 {
public:
  CellList_2(double Lx, double Ly, double x0 = 0, double y0 = 0,
    double r_cut = 1, bool comm_x = false, bool comm_y = false);

  bool has_member(int ic) const { return !cell[ic].empty(); }

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij);

  void create(std::vector<_TPar> &p_arr);
  void recreate(const std::vector<_TPar> &p_arr);
  void update(const std::vector<_TPar> &p_arr);

private:
  template <class BiFunc>
  void cell_self(BiFunc f_ij, int ic) const;

  template <class BiFunc>
  void cell_cell(BiFunc f_ij, int ic1, int ic2) const;

  std::vector<std::list<_TPar*> > cell;
  std::vector<int> list_len;
};

template <class _TPar>
class CellList_2<UniNode<_TPar>> : public _CellListBase_2 {
public:
  CellList_2(double Lx, double Ly, double x0 = 0, double y0 = 0,
    double r_cut = 1, bool comm_x = false, bool comm_y = false);
  bool has_member(int ic) const { return cell[ic] != nullptr;}

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij);

  void create(std::vector<UniNode<_TPar>> &p_arr);
  void recreate(std::vector<UniNode<_TPar>> &p_arr);
  void update(std::vector<UniNode<_TPar>> &p_arr) { recreate(p_arr); }

private:
  std::vector<UniNode<_TPar> *> cell;
};

template <class _TPar>
class CellList_2<BNode<_TPar>> : public _CellListBase_2 {
public:
  CellList_2(double Lx, double Ly, double x0 = 0, double y0 = 0,
    double r_cut = 1, bool comm_x = false, bool comm_y = false);
  bool has_member(int ic) const { return cell[ic] != nullptr; }

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij);

  void create(std::vector<BNode<_TPar>> &p_arr);
  void recreate(std::vector<BNode<_TPar>> &p_arr);
  void update(std::vector<BNode<_TPar>> &p_arr) { recreate(p_arr); }

private:
  void merge_list(int row, BNode<_TPar> **tail, BNode<_TPar> **head);

  std::vector<BNode<_TPar> *> cell;
};

template<class _TPar>
CellList_2<_TPar>::CellList_2(double Lx, double Ly, double x0, double y0,
  double r_cut, bool comm_x, bool comm_y): _CellListBase_2(Lx, Ly, x0, y0,
    r_cut, comm_x, comm_y) {
  cell.reserve(ncells);
  list_len.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.emplace_back();
    list_len.push_back(0);
  }
  std::cout << n_bins.x << "\t" << n_bins.y << "\n";
}

template <class _TPar>
CellList_2<UniNode<_TPar>>::CellList_2(double Lx, double Ly, double x0, double y0,
  double r_cut, bool comm_x, bool comm_y) : _CellListBase_2(Lx, Ly, x0, y0,
    r_cut, comm_x, comm_y) {
  cell.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.emplace_back();
  }
  std::cout << n_bins.x << "\t" << n_bins.y << "\n";
}

template <class _TPar>
CellList_2<BNode<_TPar>>::CellList_2(double Lx, double Ly, double x0, double y0,
  double r_cut, bool comm_x, bool comm_y) : _CellListBase_2(Lx, Ly, x0, y0,
    r_cut, comm_x, comm_y) {
  cell.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.emplace_back();
  }
  std::cout << n_bins.x << "\t" << n_bins.y << "\n";
}

template<class _TPar>
void CellList_2<_TPar>::create(std::vector<_TPar>& p_arr) {
  int nPar = p_arr.size();
  for (int ip = 0; ip < nPar; ip++) {
    int ic = get_ic(p_arr[ip]);
    _TPar * ptr = &p_arr[ip];
    cell[ic].push_back(ptr);
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template<class _TPar>
void CellList_2<UniNode<_TPar>>::create(std::vector<UniNode<_TPar>>& p_arr) {
  int n_node = p_arr.size();
  for (int i = 0; i < n_node; i++) {
    int ic = get_ic(p_arr[i]);
    p_arr[i].next = cell[ic];
    cell[ic] = &p_arr[i];
  }
}

template<class _TPar>
void CellList_2<BNode<_TPar>>::create(std::vector<BNode<_TPar>>& p_arr) {
  int n_node = p_arr.size();
  for (int i = 0; i < n_node; i++) {
    int ic = get_ic(p_arr[i]);
    cell[ic] = node_arr[i].append_front(cell[ic]);
  }
}

template<class _TPar>
void CellList_2<_TPar>::recreate(const std::vector<_TPar>& p_arr) {
  for (int ic = 0; ic < ncells; ic++) {
    cell[ic].clear();
  }
  create(p_arr);
}

template<class _TPar>
void CellList_2<UniNode<_TPar>>::recreate(std::vector<UniNode<_TPar>>& p_arr) {
  for (int ic = 0; ic < ncells; ic++) {
    cell[ic] = nullptr;
  }
  create(p_arr);
}

template<class _TPar>
void CellList_2<BNode<_TPar>>::recreate(std::vector<BNode<_TPar>>& p_arr) {
  for (int ic = 0; ic < ncells; ic++) {
    cell[ic] = nullptr;
  }
  create(p_arr);
}

template<class _TPar>
void CellList_2<BNode<_TPar>>::merge_list(int row, BNode<_TPar> ** tail, BNode<_TPar> ** head) {
  for (int col = 0; col < bins.x; col++) {
    if (head[col]) {
      int ic = col + row * bins.x;
      if (cell[ic]) {
        tail[col]->next = head[col];
        head[col]->prev = tail[col];
      } else {
        cell[ic] = head[col];
      }
      head[col] = nullptr;
    }
    tail[col] = nullptr;
  }
}

template<class _TPar>
void CellList_2<_TPar>::update(const std::vector<_TPar>& p_arr) {
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
template<class BiFunc>
void CellList_2<_TPar>::for_each_pair(BiFunc f_ij) {
  auto intra_cell = [this, f_ij] (int ic) {
    cell_self(f_ij, ic);
  };
  auto inter_cell = [this, f_ij](int ic1, int ic2) {
    cell_cell(f_ij, ic1, ic2);
  };
  for_each_nearest_cell(intra_cell, inter_cell, 0, n_bins.y, 0, n_bins.x);
}

template <class _TPar>
template <class BiFunc>
void CellList_2<UniNode<_TPar>>::for_each_pair(BiFunc f_ij) {
  auto intra_cell = [this, f_ij](int ic) {
    for_each_node_pair(cell[ic], f_ij);
  };
  auto inter_cell = [this, f_ij](int ic1, int ic2) {
    for_each_node_pair(cell[ic1], cell[ic2], f_ij);
  };
  for_each_nearest_cell(intra_cell, inter_cell, 0, n_bins.y, 0, n_bins.x);
}

template <class _TPar>
template <class BiFunc>
void CellList_2<BNode<_TPar>>::for_each_pair(BiFunc f_ij) {
  auto intra_cell = [this, f_ij](int ic) {
    for_each_node_pair(cell[ic], f_ij);
  };
  auto inter_cell = [this, f_ij](int ic1, int ic2) {
    for_each_node_pair(cell[ic1], cell[ic2], f_ij);
  };
  for_each_nearest_cell(intra_cell, inter_cell, 0, n_bins.y, 0, n_bins.x);
}


template<class _TPar>
template<class BiFunc>
void CellList_2<_TPar>::cell_self(BiFunc f_ij, int ic) const {
  auto end2 = cell[ic].cend();
  auto end1 = std::prev(end2);
  for (auto it1 = cell[ic].cbegin(); it1 != end1; ++it1) {
    for (auto it2 = std::next(it1); it2 != end2; ++it2) {
      f_ij(*it1, *it2);
    }
  }
}

template<class _TPar>
template<class BiFunc>
void CellList_2<_TPar>::cell_cell(BiFunc f_ij, int ic1, int ic2) const {
  auto end1 = cell[ic1].cend();
  auto end2 = cell[ic2].cend();
  auto beg2 = cell[ic2].cbegin();
  for (auto it1 = cell[ic1].cbegin(); it1 != end1; ++it1) {
    for (auto it2 = beg2; it2 != end2; ++it2) {
      f_ij(*it1, *it2);
    }
  }
}

#endif
