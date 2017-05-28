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
  void interact();
  void interact(Cell *c2);
  void interact(Cell *c2, double a, double b);
  void push_front(Node *node);
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
