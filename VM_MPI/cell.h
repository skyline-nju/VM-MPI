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
  void find_neighbor(Cell* cell, int idx, int ncols, int nrows);
  void find_neighbor(Cell* cell, int col0, int row0,
                     int ncols, int nrows, int first_row);
  void interact();
  void interact(int neighbor_idx);
  void interact(int neighbor_idx, double a, double b);
  void push_front(Node *node);
  static void update_velocity_inner_row(Cell *p, int ncols, double Lx);
  static void update_velocity_bottom_row(Cell *p, int ncols, double Lx);
  static void find_neighbor_one_row(Cell* cell, int row0,
                                    int ncols, int nrows, int first_row);
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
