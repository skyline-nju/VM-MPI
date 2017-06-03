#include "cell.h"
using namespace std;

Cell::Cell() {
  head = NULL;
  size = 0;
  for (int i = 0; i < 4; i++) {
    neighbor[i] = NULL;
  }
}

int dcol[4] = { 1, -1, 0, 1 };
int drow[4] = { 0, 1, 1, 1 };

void Cell::find_neighbor(Cell *cell, int ncols, int nrows,
                         int col0, int row0) {

  int idx = col0 + row0 * ncols;
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
    cell[idx].neighbor[i] = &cell[col + row * ncols];
  }
}

void Cell::find_all_neighbor(Cell * cell, int ncols, int nrows) {
  for (int row = 0; row < nrows; row++) {
    for (int col = 0; col < ncols; col++) {
      find_neighbor(cell, ncols, nrows, col, row);
    }
  }
}

void Cell::clear_row(int row, int ncols, Cell *cell,
                     stack<Node *> &empty_pos) {
  for (int col = 0, j = row * ncols; col < ncols; col++) {
    int idx = col + j;
    if (cell[idx].head) {
      Node *curNode = cell[idx].head;
      do {
        curNode->is_empty = true;
        empty_pos.push(curNode);
        curNode = curNode->next;
      } while (curNode);
      cell[idx].head = NULL;
      cell[idx].size = 0;
    }
  }
}

void Cell::resize(Cell **cell, int ncols, int &nrows,
                  const int *offset, stack<Node *> &empty_pos) {
  Cell *p = *cell;
  if (offset[0] == -1) {
    clear_row(0, ncols, p, empty_pos);
    nrows--;
    p += ncols;
  } else if (offset[0] == 1) {
    nrows++;
    p -= ncols;
    for (int col = 0; col < ncols; col++) {
      p[col].head = NULL;
      p[col].size = 0;
      find_neighbor(p, ncols, nrows, col, 0);
    }
  }
  if (offset[1] == -1) {
    clear_row(nrows - 1, ncols, p, empty_pos);
    nrows--;
  } else if (offset[1] == 1) {
    nrows++;
    for (int col = 0, j = (nrows - 1) * ncols; col < ncols; col++) {
      p[col + j].head = NULL;
      p[col + j].size = 0;
      find_neighbor(p, ncols, nrows, col, nrows - 2);
    }
  }
  *cell = p;
}

/****************************************************************************/
/*          Using cell list to speed up the calculation of interactions     */
/****************************************************************************/

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

void Cell::update_velocity_by_row(int row, Cell *cell, int ncols, double Lx) {
  if (row > 0)
    update_velocity_inner_row(cell + row * ncols, ncols, Lx);
  else
    update_velocity_bottom_row(cell, ncols, Lx);
}