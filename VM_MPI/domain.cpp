#include "domain.h"
#include <iostream>
using namespace std;

SubDomain::SubDomain(double _Lx, double _Ly, int tot_domain, int idx, 
                     unsigned long long seed) {
  double dLy = _Ly / tot_domain;
  Ly_l = idx * dLy - 1;
  Ly_h = (idx + 1) * dLy + 1;
  ncols = int(_Lx);
  nrows = int(dLy) + 2;
  tot_cells = ncols * nrows;
  cell = new Cell[tot_cells];
  myran = new Ran(seed + idx);
  Lx = _Lx;
  Ly = _Ly;
}

void SubDomain::update_velocity_inner() {
  for (int row = 1; row < nrows - 2; row++) {
    int j = row * ncols;
    for (int col = 0; col < ncols; col++) {
      int i = col + j;
      if (cell[i].head) {
        cell[i].interact();
        cell[i].interact(cell[i].neighbor[0]);
        cell[i].interact(cell[i].neighbor[1]);
        cell[i].interact(cell[i].neighbor[2]);
        cell[i].interact(cell[i].neighbor[3]);
      }
    }
  }
}

void SubDomain::update_velocity_edge() {
  for (int col = 0; col < ncols; col++) {
    int i = col;
    if (cell[i].head) {
      cell[i].interact(cell[i].neighbor[1]);
      cell[i].interact(cell[i].neighbor[2]);
      cell[i].interact(cell[i].neighbor[3]);
    }
  }
  for (int col = 0; col < ncols; col++) {
    int i = col + (nrows - 2) * ncols;
    if (cell[i].head) {
      cell[i].interact();
      cell[i].interact(cell[i].neighbor[0]);
      cell[i].interact(cell[i].neighbor[1]);
      cell[i].interact(cell[i].neighbor[2]);
      cell[i].interact(cell[i].neighbor[3]);
    }
  }
}

void SubDomain::update_position_inner(double eta) {
  for (auto iter = particle.begin(); iter != particle.end(); ++iter) {
    if ((*iter).valid && !(*iter).moved) {
      double noise = eta * 2 * PI * (myran->doub() - 0.5); //disorder free
      (*iter).update_coor(noise, Lx, Ly_l, Ly_h);
    }
  }
}

void SubDomain::update_position_edge(double eta) {
  // mark particles with row = 0 or nrows - 1 invalid
  for (int col = 0, j = (nrows-1) * ncols; col < ncols; col++) {
    cell[col].mark_as_invalid(empty_pos);
    cell[col + j].mark_as_invalid(empty_pos);
    cell[col + j].head = nullptr;
  }
  // row = 1
  for (int col = 0; col < ncols; col++) {
    cell[col + ncols].move(myran, eta, cell, Lx, Ly, Ly_l, Ly_h);
  }
  // row = nrows - 2
  for (int col = 0, j = (nrows - 2) * ncols; col < ncols; col++) {
    cell[col + j].move(myran, eta, cell, Lx, Ly, Ly_l, Ly_h);
  }
}

int SubDomain::get_pNum(int row) {
  int count = 0;
  for (int col = 0, j = row * ncols; col < ncols; col++) {
    count += cell[col + j].size;
  }
  return count;
}

void SubDomain::pack(int row, vector<ParticleData> &buff) {
  buff.reserve(get_pNum(row));
  for (int col = 0, j = row * ncols; col < ncols; col++) {
    int i = col + j;
    int pos = 0;
    if (cell[i].head) {
      Node *curNode = cell[i].head;
      do {
        if (!curNode->new_arrival) {
          buff.push_back(curNode->pdata);
        } else {
          curNode->new_arrival = false;
        }
        curNode = curNode->next;
      } while (curNode);
    }
  }
}

void SubDomain::unpack(const vector<ParticleData> &buff, int row) {
  int j = row * ncols;
  for (auto iter = buff.begin(); iter != buff.end(); ++iter) {
    int cell_idx = int((*iter).x) + j;
    int par_idx;
    if (!empty_pos.empty()) {
      par_idx = empty_pos.top();
      empty_pos.pop();
    } else {
      particle.push_back(Node());
      par_idx = particle.size() - 1;
    }
    particle[par_idx].pdata = *iter;
    particle[par_idx].cell_idx = cell_idx;
    particle[par_idx].par_idx = par_idx;
    particle[par_idx].valid = true;
    particle[par_idx].new_arrival = true;
    particle[par_idx].next = cell[cell_idx].head;
    cell[cell_idx].head = &particle[par_idx];
  }
}
