#ifndef NODE_H
#define NODE_H
#include "rand.h"
#include "comn.h"


struct ParticleData
{
  double x;
  double y;
  double vx;
  double vy;
};

struct Node
{
  //Node();
  //Node(const ParticleData &pdata0);
  double rr(Node *node);
  double rr(Node *node, double a, double b);
  void addV(Node *node);
  void align(Node *node);
  void align(Node *node, double a, double b);
  void move(double noise);
  void update_coor(double noise, double Lx, double Ly_l, double Ly_h);
  void update_coor(double noise, double Lx, double Ly,
                   double Ly_l, double Ly_h, bool &out_range);

  ParticleData pdata;
  double vx0;
  double vy0;
  int cell_idx;
  int par_idx;
  Node* next;
  bool moved;
  bool valid;
  bool new_arrival;

  static double Lx;
  static double Ly;
  static double v0;
  static double rho_0;
  static int N;
};

inline void Node::addV(Node *node) {
  pdata.vx += node->vx0;
  pdata.vy += node->vy0;
  node->pdata.vx += vx0;
  node->pdata.vy += vy0;
}

inline double Node::rr(Node *node) {
  double dx = node->pdata.x - pdata.x;
  double dy = node->pdata.y - pdata.y;
  return dx*dx + dy*dy;
}

inline double Node::rr(Node *node, double a, double b) {
  double dx = node->pdata.x - pdata.x + a;
  double dy = node->pdata.y - pdata.y + b;
  return dx*dx + dy*dy;
}

#endif
