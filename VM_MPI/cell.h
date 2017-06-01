#ifndef CELL_H
#define CELL_H
#include "node.h"
#include "comn.h"
#include <vector>
#include <stack>

class Cell
{
public:
  Cell();
  void interact();
  void interact(int neighbor_idx);
  void interact(int neighbor_idx, double a, double b);
  void push_front(Node *node);
  static void update_velocity_inner_row(Cell *p, int ncols, double Lx);
  static void update_velocity_bottom_row(Cell *p, int ncols, double Lx);
  static void find_neighbor(Cell* cell, int ncols, int nrows,
                            int mycol, int myrow);
  static void find_all_neighbor(Cell *cell, int ncols, int nrows);
  static void update_neighbor(Cell **cell, int ncols, int &nrows,
                              const int *offset);
  Node* head;
  int size;
  Cell *neighbor[4];
};

inline void Cell::push_front(Node *node) {
  node->next = head;
  head = node;
  size++;
}

#endif
