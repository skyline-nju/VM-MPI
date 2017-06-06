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
  int cell_idx(double yl, int ncols);
  void update_coor(double noise, double Lx, double yl);
  static void ini_random(Node *par, int nPar, Ran *myran,
                         double Lx, double yl, double yh);
  static void ini_from_snap(Node **par, int &npar, int &max_par_num,
                            double magnification, const std::string &filename,
                            double Lx, double Ly, int tot_rank, int myrank);
  static void sum_v(const Node *par, int end_pos,
                    int &npar, double &svx, double &svy);
  static int get_nPar(const Node *par, int end_pos);

  double x;
  double y;
  double vx;
  double vy;
  double vx0;
  double vy0;
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

inline int Node::cell_idx(double yl, int ncols) {
  return int(x) + int(y - yl) * ncols;
}
#endif
