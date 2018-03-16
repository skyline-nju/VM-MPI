#include "cellList.h"

#ifndef USE_MPI
CellListBase_2::CellListBase_2(double Lx, double Ly, double rcut) {
  bins.x = int(Lx / rcut);
  bins.y = int(Ly / rcut);
  l_cell.x = Lx / bins.x;
  l_cell.y = Ly / bins.y;
  ncells = bins.x * bins.y;
}

CellListIdx_2::CellListIdx_2(double Lx, double Ly, double rcut)
  : CellListBase_2(Lx, Ly, rcut) {
  cell.reserve(ncells);
  list_len.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.emplace_back();
    list_len.push_back(0);
  }
}

CellListNode_2::CellListNode_2(double Lx, double Ly, double rcut)
  : CellListBase_2(Lx, Ly, rcut) {
  cell.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.push_back(nullptr);
  }
}

CellListBiNode_2::CellListBiNode_2(double Lx, double Ly, double rcut)
  : CellListBase_2(Lx, Ly, rcut) {
  std::cout << "ncells = " << ncells << std::endl;
  cell.reserve(ncells);
  //tail.reserve(ncells);
  //head2.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.push_back(nullptr);
    //tail.push_back(nullptr);
    //head2.push_back(nullptr);
  }
}

#else
CellListBase_2::CellListBase_2(double Lx, double Ly, double y0, double r_cut) {
  bins.x = int(Lx / r_cut);
  bins.y = int(Ly / r_cut) + 2;
  l_cell.x = Lx / bins.x;
  l_cell.y = Ly / (bins.y - 2);
  ncells = bins.x * bins.y;
  y_l = y0;
  y_h = y_l + Ly;
  y_l_ext = y_l - l_cell.y;
  y_h_ext = y_h + l_cell.y;
}

CellListIdx_2::CellListIdx_2(double Lx, double Ly, double y0, double r_cut)
  : CellListBase_2(Lx, Ly, y0, r_cut) {
  cell.reserve(ncells);
  list_len.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.emplace_back();
    list_len.push_back(0);
  }
}

CellListNode_2::CellListNode_2(double Lx, double Ly, double y0, double r_cut)
  : CellListNode_2(Lx, Ly, y0, r_cut) {
  cell.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.push_back(nullptr);
  }
}

CellListBiNode_2::CellListBiNode_2(double Lx, double Ly, double y0, double rcut)
  : CellListBase_2(Lx, Ly, y0, rcut) {
  head.reserve(ncells);
  tail.reserve(ncells);
  length.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    head.push_back(nullptr);
    tail.push_back(nullptr);
    length.push_back(0);
  }
}

#endif

void CellListNode_2::create(std::vector<Node> &node_arr) {
  int n_node = node_arr.size();
  for (int i = 0; i < n_node; i++) {
    int ic = get_ic(node_arr[i]);
    node_arr[i].next = cell[ic];
    cell[ic] = &node_arr[i];
  }
}

void CellListNode_2::recreate(std::vector<Node>& node_arr) {
  for (int ic = 0; ic < ncells; ic++) {
    cell[ic] = nullptr;
  }
  create(node_arr);
}

void CellListBiNode_2::create(std::vector<BiNode> &node_arr) {
  int n_node = node_arr.size();
  for (int i = 0; i < n_node; i++) {
    int ic = get_ic(node_arr[i]);
    cell[ic] = node_arr[i].append_front(cell[ic]);
  }
}

void CellListBiNode_2::recreate(std::vector<BiNode>& node_arr) {
  for (int ic = 0; ic < ncells; ic++) {
    cell[ic] = nullptr;
  }
  create(node_arr);
}

void CellListBiNode_2::merge_list(int row, BiNode ** tail, BiNode ** head) {
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

