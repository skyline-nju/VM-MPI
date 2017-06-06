/****************************************************************************/
/*          Using cell list to speed up the calculation of interactions     */
/*          between two particles within range 1.                           */
/****************************************************************************/
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
  static void update_velocity_by_row(int row, 
                                     Cell *cell, int ncols, double Lx);
  static void find_neighbor(Cell* cell, int ncols, int nrows,
                            int mycol, int myrow);
  static void find_all_neighbor(Cell *cell, int ncols, int nrows);
  static void clear_row(int row, int ncols, Cell *cell,
                        std::stack<Node *> &empty_pos);
  static void resize(Cell **cell, int ncols, int &nrows,
                     const int *offset, std::stack<Node *> &empty_pos);
  static int get_nPar(const Cell *cell, int ncols, int nrows);
  static void get_nPar(const Cell *cell, int ncols, int nrows, int *count);

 
  Node* head;
  int size;
  double disorder;
  Cell *neighbor[4];    // right, upper left, upper, upper right
};

inline void Cell::push_front(Node *node) {
  node->next = head;
  head = node;
  size++;
}

inline void Cell::update_velocity_by_row(int row, Cell *cell, int ncols, double Lx) {
  if (row > 0)
    update_velocity_inner_row(cell + row * ncols, ncols, Lx);
  else
    update_velocity_bottom_row(cell, ncols, Lx);
}

#endif