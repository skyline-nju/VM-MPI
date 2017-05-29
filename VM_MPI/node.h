#ifndef NODE_H
#define NODE_H
#include "rand.h"
#include "comn.h"

struct Node
{
  double rr(Node *node);
  double rr(Node *node, double a, double b);
  void addV(Node *node);
  void align(Node *node);
  void align(Node *node, double a, double b);
  void cal_cell_idx(double yl, int ncols);
  void update_coor(double noise, double Lx, double yl);

  double x;
  double y;
  double vx;
  double vy;
  double vx0;
  double vy0;
  int cell_idx;
  int par_idx;
  Node* next;
  bool is_moved;
  bool is_empty;
  bool is_ghost;
  bool new_arrival;

  static double v0;
};

inline void Node::addV(Node *node) {
  vx += node->vx0;
  vy += node->vy0;
  node->vx += vx0;
  node->vy += vy0;
}

inline double Node::rr(Node *node) {
  double dx = node->x - x;
  double dy = node->y - y;
  return dx*dx + dy*dy;
}

inline double Node::rr(Node *node, double a, double b) {
  double dx = node->x - x + a;
  double dy = node->y - y + b;
  return dx*dx + dy*dy;
}

inline void Node::cal_cell_idx(double yl, int ncols) {
  int col = int(x);
  int row = int(y - yl);
  cell_idx = col + row * ncols;
}
#endif
