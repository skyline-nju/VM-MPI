#ifndef CELL_H
#define CELL_H
#include "node.h"
#include "comn.h"
#include <vector>
#include <stack>

class Cell
{
public:
  Cell() { head = nullptr;  size = 0; }
  void interact();
  void interact(Cell *c2);
  void mark_as_invalid(std::stack<unsigned int> &empty_pos);
  void move(Ran *myran, double eta, Cell *cell, double Lx, double Ly,
            double Ly_l, double Ly_h);

  Node* head;
  int size;
  std::vector<Cell *> neighbor;
};
#endif
