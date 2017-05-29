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

void Cell::interact(Cell *c2) {
  if (c2->head) {
    Node *node1 = head;
    Node *node2;
    do {
      node2 = c2->head;
      do {
        node1->align(node2);
        node2 = node2->next;
      } while (node2);
      node1 = node1->next;
    } while (node1);
  }
}

void Cell::interact(Cell *c2, double a, double b) {
  if (c2->head) {
    Node *node1 = head;
    Node *node2;
    do {
      node2 = c2->head;
      do {
        node1->align(node2, a, b);
        node2 = node2->next;
      } while (node2);
      node1 = node1->next;
    } while (node1);
  }
}

