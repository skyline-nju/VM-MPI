#ifndef GRID_H
#define GRID_H
#include "node.h"

struct Grid
{
  Grid() { head = nullptr; }
  void cell_cell();
  void cell_cell(Grid*);
  void cell_cell(Grid*, double, double);
  static void all_pairs(Grid *cell);
  static void link_nodes(Grid *cell, Node *node);
  static void refresh(Grid *cell, Node *node);
  static Grid *ini(double Lx, double Ly);

  Node* head;

  static int mx;
  static int my;
  static int mm;
  static double l0;
};
#endif


