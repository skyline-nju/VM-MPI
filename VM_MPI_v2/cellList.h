#ifndef CELLLIST_H
#define CELLLIST_H
#include <vector>
#include <list>
#include <iostream>
#include "vect.h"
#include "comn.h"
#include "particle.h"

class CellListBase_2 {
public:
#ifndef USE_MPI
  CellListBase_2(double Lx, double Ly, double rcut);
#else
  CellListBase_2(double Lx, double Ly, double y0, double rcut = 1);
#endif

  template <class Par>
  int get_ic(const Par &p) const;

  virtual bool has_member(int ic) const = 0;

  template <class UniFunc>
  void for_each_nearest_cell(UniFunc f, int row) const;
protected:
  Vec_2<double> l_cell;
  Vec_2<int> bins;
  int ncells;
#ifdef USE_MPI
  double y_l;
  double y_h;
  double y_l_ext;
  double y_h_ext;
#endif
};

template <class Par>
int CellListBase_2::get_ic(const Par &p) const {
#ifndef USE_MPI
  return int(p.x / l_cell.x) + int(p.y / l_cell.y) * bins.x;
#else
  return int(p.x / l_cell.x) + int((p.y - y_l_ext) / l_cell.y) * bins.x;
#endif
}

template<class UniFunc>
void CellListBase_2::for_each_nearest_cell(UniFunc f1, int row) const {
  int idx_cell[5];
  int iy_up = row + 1;
  if (iy_up == bins.y)
    iy_up = 0;
  int iy_times_nx = row * bins.x;
  int iy_up_times_nx = iy_up * bins.x;
  for (int ix = 0; ix < bins.x; ix++) {
    int my_i = ix + iy_times_nx;
    if (has_member(my_i)) {
      int ix_rt = ix + 1;
      if (ix_rt == bins.x)
        ix_rt = 0;
      int ix_lt = ix - 1;
      if (ix_lt == -1)
        ix_lt = bins.x - 1;
      idx_cell[0] = my_i;
      idx_cell[1] = ix_rt + iy_times_nx;
      idx_cell[2] = ix_lt + iy_up_times_nx;
      idx_cell[3] = ix + iy_up_times_nx;
      idx_cell[4] = ix_rt + iy_up_times_nx;
      f1(idx_cell);
    }
  }
}

class CellListIdx_2 : public CellListBase_2 {
public:
#ifndef USE_MPI
  CellListIdx_2(double Lx, double Ly, double rcut);
#else
  CellListIdx_2(double Lx, double Ly, double y0, double rcut = 1);
#endif
  bool has_member(int ic) const { return !cell[ic].empty(); }

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij, int mode = 0) const;

  template <class BiFunc>
  void cell_cell(BiFunc f_ij, int ic) const;

  template <class BiFunc>
  void cell_cell(BiFunc f_ij, int ic1, int ic2) const;

  template <class Par>
  void create(const std::vector<Par> &p_arr);

  template <class Par>
  void recreate(const std::vector<Par> &p_arr);

  template <class Par>
  void update(const std::vector<Par> &p_arr);

  template <class Par, class BiFunc>
  void cal_force(const std::vector<Par> &p_arr, BiFunc f_ij);
#ifdef BINARY_PARTICLE

  template <class Par1, class Par2>
  void create(const std::vector<Par1> &p1_arr,
    const std::vector<Par2> &p2_arr);
  template <class Par1, class Par2>
  void recreate(const std::vector<Par1> &p1_arr,
    const std::vector<Par2> &p2_arr);

  template <class Par1, class Par2>
  void update(const std::vector<Par1> &p1_arr,
    const std::vector<Par2> &p2_arr);
  template <class Par1, class Par2, class BiFunc>
  void cal_force(const std::vector<Par1> &p1_arr,
    const std::vector<Par2> &p2_arr, BiFunc f_ij);
#endif

#ifdef SPATIAL_SORT
  template <class Par, class BiFunc>
  void cal_force(std::vector<Par>, BiFunc pair_force,
    MySpatialSortingTraits<Par> &sst);
#endif


protected:
  std::vector< std::list<int> > cell;
  std::vector<int> list_len;
  int count;
};

template<class BiFunc>
void CellListIdx_2::for_each_pair(BiFunc f_ij, int mode) const {
  auto lambda = [this, f_ij](int *ic) {
    if (cell[ic[0]].size() > 1)
      cell_cell(f_ij, ic[0]);
    if (!cell[ic[1]].empty())
      cell_cell(f_ij, ic[0], ic[1]);
    if (!cell[ic[2]].empty())
      cell_cell(f_ij, ic[0], ic[2]);
    if (!cell[ic[3]].empty())
      cell_cell(f_ij, ic[0], ic[3]);
    if (!cell[ic[4]].empty())
      cell_cell(f_ij, ic[0], ic[4]);
  };
  if (mode == 0) {
    for (int row = 0; row < bins.y; row++) {
      for_each_nearest_cell(lambda, row);
    }
  } else if (mode == 1) {
    for (int row = 1; row < bins.y - 2; row++) {
      for_each_nearest_cell(lambda, row);
    }
  } else {
    for_each_nearest_cell(lambda, bins.y - 2);
    for_each_nearest_cell(lambda, 0);
  }
}

template<class BiFunc>
void CellListIdx_2::cell_cell(BiFunc f_ij, int ic) const {
  auto end2 = cell[ic].cend();
  auto end1 = std::prev(end2);
  for (auto it1 = cell[ic].cbegin(); it1 != end1; ++it1) {
    for (auto it2 = std::next(it1); it2 != end2; ++it2) {
      f_ij(*it1, *it2);
    }
  }
}

template<class BiFunc>
void CellListIdx_2::cell_cell(BiFunc f_ij, int ic1, int ic2) const {
  auto end1 = cell[ic1].cend();
  auto end2 = cell[ic2].cend();
  auto beg2 = cell[ic2].cbegin();
  for (auto it1 = cell[ic1].cbegin(); it1 != end1; ++it1) {
    for (auto it2 = beg2; it2 != end2; ++it2) {
      f_ij(*it1, *it2);
    }
  }
}

template<class Par>
void CellListIdx_2::create(const std::vector<Par>& p_arr) {
  int nPar = p_arr.size();
  for (int ip = 0; ip < nPar; ip++) {
    cell[get_ic(p_arr[ip])].push_back(ip);
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template<class Par>
void CellListIdx_2::recreate(const std::vector<Par>& p_arr) {
  for (int ic = 0; ic < ncells; ic++) {
    cell[ic].clear();
  }
  create(p_arr);
}

template<class Par>
void CellListIdx_2::update(const std::vector<Par>& p_arr) {
  for (int iy = 0; iy < bins.y; iy++) {
    int iy_times_nx = iy * bins.x;

    for (int ix = 0; ix < bins.x; ix++) {
      int ic = ix + iy_times_nx;
      int depth = list_len[ic];
      if (depth) {
        int it_count = 0;
        for (auto it = cell[ic].begin(); it_count < depth; it_count++) {
          int ip = *it;
          int ic_new = get_ic(p_arr[ip]);
          if (ic_new == ic) {
            ++it;
          } else {
            it = cell[ic].erase(it);
            //if (ic_new >= cell.size()) {
            //  std::cout << ic_new << std::endl;
            //}
            cell[ic_new].push_back(ip);
          }
        }
      }
    }
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template<class Par, class BiFunc>
inline void CellListIdx_2::cal_force(const std::vector<Par>& p_arr, BiFunc f_ij) {
  update(p_arr);
  for_each_pair(f_ij);
}

#ifdef BINARY_PAR
template<class Par1, class Par2>
void CellListIdx_2::create(const std::vector<Par1>& p1_arr, const std::vector<Par2> &p2_arr) {
  int n1 = p1_arr.size();
  for (int ip = 0; ip < n1; ip++) {
    cell[get_ic(p1_arr[ip])].push_back(ip);
  }
  int n2 = p2_arr.size();
  for (int ip = 0; ip < n2; ip++) {
    cell[get_ic(p2_arr[ip])].push_back(ip + n1);
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template <class Par1, class Par2>
void CellListIdx_2::recreate(const std::vector<Par1> &p1_arr, const std::vector<Par2> &p2_arr) {
  for (int ic = 0; ic < ncells; ic++) {
    cell[ic].clear();
  }
  create(p1_arr, p2_arr);
}

template <class Par1, class Par2>
void CellListIdx_2::update(const std::vector<Par1> &p1_arr, const std::vector<Par2> &p2_arr) {
  int n1 = p1_arr.size();
  for (int iy = 0; iy < bins.y; iy++) {
    int iy_times_nx = iy * bins.x;

    for (int ix = 0; ix < bins.x; ix++) {
      int ic = ix + iy_times_nx;
      int depth = list_len[ic];
      if (depth) {
        int it_count = 0;
        for (auto it = cell[ic].begin(); it_count < depth; it_count++) {
          int ip = *it;
          int ic_new = ip < n1 ? get_ic(p1_arr[ip]) : get_ic(p2_arr[ip - n1]);
          if (ic_new == ic) {
            ++it;
          } else {
            it = cell[ic].erase(it);
            cell[ic_new].push_back(ip);
          }
        }
      }
    }
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template <class Par1, class Par2, class BiFunc>
inline void CellListIdx_2::cal_force(const std::vector<Par1> &p1_arr,
                                     const std::vector<Par2> &p2_arr, BiFunc f_ij) {
  update(p1_arr, p2_arr);
  for_each_pair(f_ij);
}
#endif

#ifdef SPATIAL_SORT
template<class Par, class BiFunc>
inline void CellListIdx_2::cal_force(std::vector<Par> p_arr, BiFunc pair_force,
  MySpatialSortingTraits<Par>& sst) {
  if (count % 1000 == 0) {
    CGAL::spatial_sort(p_arr.begin(), p_arr.end(), sst);
    recreate(p_arr);
  } else {
    update(p_arr);
  }
  for_each_pair(pair_force);
}
#endif

class CellListNode_2 : public CellListBase_2 {
public:
#ifndef USE_MPI
  CellListNode_2(double Lx, double Ly, double rcut);
#else
  CellListNode_2(double Lx, double Ly, double y0, double rcut = 1);
#endif

  bool has_member(int ic) const { return cell[ic]; }

  void create(std::vector<Node> &node_arr);

  void recreate(std::vector<Node> &node_arr);

  template <class BiFunc>
  void cell_cell(BiFunc f_ij, int ic) const;

  template <class BiFunc>
  void cell_cell(BiFunc f_ij, int ic1, int ic2) const;

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij, int mode = 0) const;

  template <class UniFunc>
  void for_each(UniFunc f_i) const;

protected:
  std::vector<Node *> cell;
};

template<class BiFunc>
void CellListNode_2::cell_cell(BiFunc f_ij, int ic) const {
  Node *node1 = cell[ic];
  Node *node2;
  while (node1->next) {
    node2 = node1->next;
    do {
      f_ij(node1, node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  }
}

template<class BiFunc>
void CellListNode_2::cell_cell(BiFunc f_ij, int ic1, int ic2) const {
  Node* node1 = cell[ic1];
  Node* node2;
  do {
    node2 = cell[ic2];
    do {
      f_ij(node1, node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  } while (node1);
}
template<class BiFunc>
void CellListNode_2::for_each_pair(BiFunc f_ij, int mode) const {
  auto lambda = [this, f_ij](int *ic) {
    cell_cell(f_ij, ic[0]);
    if (cell[ic[1]])
      cell_cell(f_ij, ic[0], ic[1]);
    if (cell[ic[2]])
      cell_cell(f_ij, ic[0], ic[2]);
    if (cell[ic[3]])
      cell_cell(f_ij, ic[0], ic[3]);
    if (cell[ic[4]])
      cell_cell(f_ij, ic[0], ic[4]);
  };
  if (mode == 0) {
    for (int row = 0; row < bins.y; row++) {
      for_each_nearest_cell(lambda, row);
    }
  } else if (mode == 1) {
    for (int row = 1; row < bins.y - 2; row++) {
      for_each_nearest_cell(lambda, row);
    }
  } else {
    for_each_nearest_cell(lambda, bins.y - 2);
    for_each_nearest_cell(lambda, 0);
  }
}

template<class UniFunc>
void CellListNode_2::for_each(UniFunc f_i) const {
  for (int ic = 0; ic < ncells; ic++) {
    if (cell[ic]) {
      Node * cur = cell[ic];
      do {
        f_i(cur);
        cur = cur->next;
      } while (cur);
    }
  }
}


class CellListBiNode_2 : public CellListBase_2 {
public:
#ifndef USE_MPI
  CellListBiNode_2(double Lx, double Ly, double rcut);
#else
  CellListBiNode_2(double Lx, double Ly, double y0, double rcut = 1);
#endif

  bool has_member(int ic) const { return cell[ic]; }

  void create(std::vector<BiNode> &node_arr);

  void recreate(std::vector<BiNode> &node_arr);

  template <class BiFunc>
  void cell_cell(BiFunc f_ij, int ic) const;

  template <class BiFunc>
  void cell_cell(BiFunc f_ij, int ic1, int ic2) const;

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij, int mode = 0) const;

  template <class UniFunc>
  void for_each(UniFunc f_i) const;

  template <class UniFunc>
  void update(UniFunc f_move);

  void merge_list(int row, BiNode ** tail, BiNode ** head);

  template <class BiFunc>
  void update_by_row(BiFunc f_move);

  template <class BiFunc>
  void update(BiFunc f_move, int row, BiNode **tail_cur, BiNode **head_pre,
    BiNode **head_cur, BiNode **head_next);

  void leave_cell(BiNode** head, BiNode **curNode, int ic, int new_col);


protected:
  std::vector<BiNode *> cell;
};

template<class BiFunc>
void CellListBiNode_2::cell_cell(BiFunc f_ij, int ic) const {
  BiNode *node1 = cell[ic];
  BiNode *node2;
  while (node1->next) {
    node2 = node1->next;
    do {
      f_ij(node1, node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  }
}

template<class BiFunc>
void CellListBiNode_2::cell_cell(BiFunc f_ij, int ic1, int ic2) const {
  BiNode* node1 = cell[ic1];
  BiNode* node2;
  do {
    node2 = cell[ic2];
    do {
      f_ij(node1, node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  } while (node1);
}

template<class BiFunc>
void CellListBiNode_2::for_each_pair(BiFunc f_ij, int mode) const {
  auto lambda = [this, f_ij](int *ic) {
    cell_cell(f_ij, ic[0]);
    if (cell[ic[1]])
      cell_cell(f_ij, ic[0], ic[1]);
    if (cell[ic[2]])
      cell_cell(f_ij, ic[0], ic[2]);
    if (cell[ic[3]])
      cell_cell(f_ij, ic[0], ic[3]);
    if (cell[ic[4]])
      cell_cell(f_ij, ic[0], ic[4]);
  };
  if (mode == 0) {
    for (int row = 0; row < bins.y; row++) {
      for_each_nearest_cell(lambda, row);
    }
  } else if (mode == 1) {
    for (int row = 1; row < bins.y - 2; row++) {
      for_each_nearest_cell(lambda, row);
    }
  } else {
    for_each_nearest_cell(lambda, bins.y - 2);
    for_each_nearest_cell(lambda, 0);
  }
}

template<class UniFunc>
void CellListBiNode_2::for_each(UniFunc f_i) const {
  for (int ic = 0; ic < ncells; ic++) {
    if (cell[ic]) {
      BiNode * cur = cell[ic];
      do {
        f_i(cur);
        cur = cur->next;
      } while (cur);
    }
  }
}

template<class UniFunc>
void CellListBiNode_2::update(UniFunc f_move) {
  std::vector<BiNode*> tail(ncells, nullptr);
  std::vector<BiNode*> head(ncells, nullptr);
  for (int ic = 0; ic < ncells; ic++) {
    if (cell[ic]) {
      BiNode * cur_node = cell[ic];
      do {
        f_move(cur_node);
        int ic_new = get_ic(*cur_node);
        if (ic_new != ic) {
          cur_node->break_away(&cell[ic]);
          BiNode* new_node = cur_node;
          cur_node = cur_node->next;
          new_node->append_front(&head[ic_new]);
        } else {
          tail[ic] = cur_node;
          cur_node = cur_node->next;
        }
      } while (cur_node);
    }
  }
  for (int ic = 0; ic < ncells; ic++) {
    if (head[ic]) {
      if (cell[ic]) {
        tail[ic]->next = head[ic];
        head[ic]->prev = tail[ic];
      } else {
        cell[ic] = head[ic];
      }
    }
  }
}

template<class BiFunc>
void CellListBiNode_2::update_by_row(BiFunc f_move) {
  BiNode ** tail_0 = new BiNode*[bins.x]{};
  BiNode ** tail_1 = new BiNode*[bins.x]{};
  BiNode ** tail_2 = new BiNode*[bins.x]{};
  BiNode ** head_0 = new BiNode*[bins.x]{};
  BiNode ** head_1 = new BiNode*[bins.x]{};
  BiNode ** head_2 = new BiNode*[bins.x]{};
  BiNode ** head_3 = new BiNode*[bins.x]{};
  BiNode ** head_4 = new BiNode*[bins.x]{};

  auto swap2 = [](BiNode *** a, BiNode *** b) {
    BiNode ** tmp = *b;
    *b = *a;
    *a = tmp;
  };

  update(f_move, 0, tail_0, head_4, head_0, head_1);
  update(f_move, 1, tail_1, head_0, head_1, head_2);

  for (int row = 2; row < bins.y - 2; row++) {
    update(f_move, row, tail_2, head_1, head_2, head_3);
    merge_list(row - 1, tail_1, head_1);
    swap2(&tail_1, &tail_2);
    swap2(&head_1, &head_3);
    swap2(&head_1, &head_2);
  }

  update(f_move, bins.y - 2, tail_2, head_1, head_2, head_4);
  merge_list(bins.y - 3, tail_1, head_1);
  swap2(&tail_1, &tail_2);

  update(f_move, bins.y - 1, tail_2, head_2, head_4, head_0);
  merge_list(bins.y - 2, tail_1, head_2);

  merge_list(bins.y - 1, tail_2, head_4);
  merge_list(0, tail_0, head_0);

  delete[] tail_0;
  delete[] tail_1;
  delete[] tail_2;
  delete[] head_0;
  delete[] head_1;
  delete[] head_2;
  delete[] head_3;
  delete[] head_4;
}

template<class BiFunc>
void CellListBiNode_2::update(BiFunc f_move, int row, BiNode ** tail_cur,
                              BiNode ** head_pre, BiNode ** head_cur, BiNode ** head_next) {
  for (int col = 0; col < bins.x; col++) {
    int ic = col + row * bins.x;
    if (cell[ic]) {
      Vec_2<int> d_cell;
      BiNode * cur_node = cell[ic];
      do {
        f_move(cur_node, d_cell);
        if (d_cell.y == -1) {
          leave_cell(head_pre, &cur_node, ic, col + d_cell.x);
        } else if (d_cell.y == 1) {
          leave_cell(head_next, &cur_node, ic, col + d_cell.x);
        } else if (d_cell.x != 0) {
          leave_cell(head_cur, &cur_node, ic, col + d_cell.x);
        } else {
          tail_cur[col] = cur_node;
          cur_node = cur_node->next;
        }
      } while (cur_node);
    }
  }
}

inline void CellListBiNode_2::leave_cell(BiNode** head, BiNode **curNode, int ic, int new_col) {
  (*curNode)->break_away(&cell[ic]);
  BiNode * new_node = *curNode;
  *curNode = (*curNode)->next;
  if (new_col < 0)
    new_col += bins.x;
  else if (new_col >= bins.x)
    new_col = 0;
  new_node->append_front(&head[new_col]);
}

#endif
