#include "cell.h"
using namespace std;

Cell::Cell() {
  head = NULL;
  size = 0;
  for (int i = 0; i < 4; i++) {
    neighbor[i] = NULL;
  }
}

void Cell::find_neighbor(Cell *cell, int idx, int ncols, int nrows) {
  int col0 = idx % ncols;
  int row0 = idx / ncols;
  int dcol[4] = { 1, -1, 0, 1 };
  int drow[4] = { 0, 1, 1, 1 };
  for (int i = 0; i < 4; i++) {
    int col = col0 + dcol[i];
    int row = row0 + drow[i];
    if (col < 0) {
      col += ncols;
    } else if (col >= ncols) {
      col -= ncols;
    }
    if (row >= nrows) {
      row -= nrows;
    }
    int idx_neighbor = col + row * ncols;
    neighbor[i] = &cell[idx_neighbor];
  }
}

void Cell::find_neighbor(Cell *cell, int col0, int row0,
                         int ncols, int nrows, int first_row) {
  int dcol[4] = { 1, -1, 0, 1 };
  int drow[4] = { 0, 1, 1, 1 };
  for (int i = 0; i < 4; i++) {
    int col = col0 + dcol[i];
    int row = row0 + drow[i];
    if (col < 0) {
      col += ncols;
    } else if (col >= ncols) {
      col -= ncols;
    }
    if (row < first_row) {
      row += nrows;
    } else if (row >= first_row + nrows) {
      row -= nrows;
    }
    int idx_neighbor = col + row * ncols;
    neighbor[i] = cell + idx_neighbor;
  }
}

void Cell::interact() {
  Node *node1 = head;
  Node *node2;
  while (node1->next) {
    node2 = node1->next;
    do {
      node1->align(node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  }
}

void Cell::interact(int neighbor_idx) {
  if (neighbor[neighbor_idx]->head) {
    Node *node1 = head;
    Node *node2;
    do {
      node2 = neighbor[neighbor_idx]->head;
      do {
        node1->align(node2);
        node2 = node2->next;
      } while (node2);
      node1 = node1->next;
    } while (node1);
  }
}

void Cell::interact(int neighbor_idx, double a, double b) {
  if (neighbor[neighbor_idx]->head) {
    Node *node1 = head;
    Node *node2;
    do {
      node2 = neighbor[neighbor_idx]->head;
      do {
        node1->align(node2, a, b);
        node2 = node2->next;
      } while (node2);
      node1 = node1->next;
    } while (node1);
  }
}

void Cell::update_velocity_inner_row(Cell *p, int ncols, double Lx) {
  /* leftmost column 0 */
  {
    if (p->head) {
      p->interact();
      p->interact(0);
      p->interact(1, -Lx, 0);
      p->interact(2);
      p->interact(3);
    }
    p++;
  }
  /* inner colmuns from 1 to ncols - 2 */
  for (int col = 1; col < ncols - 1; col++) {
    if (p->head) {
      p->interact();
      p->interact(0);
      p->interact(1);
      p->interact(2);
      p->interact(3);
    }
    p++;
  }
  /* rightmost column ncols - 1 */
  {
    if (p->head) {
      p->interact();
      p->interact(0, Lx, 0);
      p->interact(1);
      p->interact(2);
      p->interact(3, Lx, 0);
    }
  }
}

void Cell::update_velocity_bottom_row(Cell *p, int ncols, double Lx) {
  /* leftmost column 0 */
  {
    if (p->head) {
      p->interact(1, -Lx, 0);
      p->interact(2);
      p->interact(3);
    }
    p++;
  }
  /* inner colmuns from 1 to ncols - 2 */
  for (int col = 1; col < ncols - 1; col++) {
    if (p->head) {
      p->interact(1);
      p->interact(2);
      p->interact(3);
    }
  }
  /* rightmost column ncols - 1 */
  {
    int i = ncols - 1;
    if (p->head) {
      p->interact(1);
      p->interact(2);
      p->interact(3, Lx, 0);
    }
  }
}

void Cell::find_neighbor_one_row(Cell * cell, int row0,
                                 int ncols, int nrows, int first_row) {
  for (int col = 0; col < ncols; col++) {
    cell[col + row0 * ncols].find_neighbor(
      cell, col, row0, ncols, nrows, first_row);
  }
}

