#include "node.h"

using namespace std;
double Node::v0 = 0.5;

void Node::align(Node *node) {
  if (rr(node) < 1)
    addV(node);
}

void Node::align(Node *node, double a, double b) {
  if (rr(node, a, b) < 1)
    addV(node);
}

void Node::update_coor(double noise, double Lx, double yl) {
  double tmp = sqrt(vx*vx + vy*vy);
  double c1 = vx / tmp;
  double s1 = vy / tmp;
  double c2 = cos(noise);
  double s2 = sin(noise);
  vx = vx0 = c1 * c2 - s1 * s2;
  vy = vy0 = c1 * s2 + c2 * s1;
  x += v0*vx;
  if (x >= Lx) {
    x -= Lx;
  } else if (x < 0) {
    x += Lx;
  }
  y += v0*vy;
  is_moved = true;
  new_arrival = false;
  cal_cell_idx(yl, int(Lx));
}

