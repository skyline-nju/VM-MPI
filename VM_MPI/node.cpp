#include "node.h"

using namespace std;
double Node::Lx;
double Node::Ly;
double Node::v0 = 0.5;
double Node::rho_0;
int Node::N;

//Node::Node() {
//  pdata.x = pdata.y = pdata.vx = pdata.vy = 0;
//  vx0 = vy0 = 0;
//  cell_idx = par_idx = 0;
//  next = nullptr;
//  moved = false;
//  valid = true;
//  new_arrival = false;
//}
//
//Node::Node(const ParticleData &pdata0) {
//  pdata = pdata0;
//  vx0 = vy0 = 0;
//  cell_idx = par_idx = 0;
//  next = nullptr;
//  moved = false;
//  valid = true;
//  new_arrival = true;
//}

void Node::align(Node *node) {
  if (rr(node) < 1)
    addV(node);
}

void Node::align(Node *node, double a, double b) {
  if (rr(node, a, b) < 1)
    addV(node);
}

void Node::move(double noise) {
  double tmp = sqrt(pdata.vx*pdata.vx + pdata.vy*pdata.vy);
  double c1 = pdata.vx / tmp;
  double s1 = pdata.vy / tmp;
  double c2 = cos(noise);
  double s2 = sin(noise);
  pdata.vx = vx0 = c1 * c2 - s1 * s2;
  pdata.vy = vy0 = c1 * s2 + c2 * s1;

  pdata.x += v0*pdata.vx;
  if (pdata.x >= Lx)
    pdata.x -= Lx;
  else if (pdata.x < 0)
    pdata.x += Lx;
  pdata.y += v0*pdata.vy;
  if (pdata.y >= Ly)
    pdata.y -= Ly;
  else if (pdata.y < 0)
    pdata.y += Ly;

}

void Node::update_coor(double noise, double Lx, double Ly_l, double Ly_h) {
  double tmp = sqrt(pdata.vx*pdata.vx + pdata.vy*pdata.vy);
  double c1 = pdata.vx / tmp;
  double s1 = pdata.vy / tmp;
  double c2 = cos(noise);
  double s2 = sin(noise);
  pdata.vx = vx0 = c1 * c2 - s1 * s2;
  pdata.vy = vy0 = c1 * s2 + c2 * s1;
  pdata.x += v0*pdata.vx;
  if (pdata.x >= Lx) {
    pdata.x -= Lx;
  } else if (pdata.x < 0) {
    pdata.x += Lx;
  }
  pdata.y += v0*pdata.vy;
  int col = int(pdata.x);
  int row = int(pdata.y - Ly_l);
  cell_idx = col + row * int(Lx);
  moved = true;
}

void Node::update_coor(double noise, double Lx, double Ly,
                       double Ly_l, double Ly_h, bool &out_range) {
  double tmp = sqrt(pdata.vx*pdata.vx + pdata.vy*pdata.vy);
  double c1 = pdata.vx / tmp;
  double s1 = pdata.vy / tmp;
  double c2 = cos(noise);
  double s2 = sin(noise);
  pdata.vx = vx0 = c1 * c2 - s1 * s2;
  pdata.vy = vy0 = c1 * s2 + c2 * s1;
  pdata.x += v0*pdata.vx;
  if (pdata.x >= Lx) {
    pdata.x -= Lx;
  } else if (pdata.x < 0) {
    pdata.x += Lx;
  }
  pdata.y += v0*pdata.vy;
  if (pdata.y < Ly_l + 1 || pdata.y >= Ly_h - 1) {
    out_range = true;
  } else {
    out_range = false;
  }
  int col = int(pdata.x);
  int row = int(pdata.y - Ly_l);
  cell_idx = col + row * int(Lx);
  if (pdata.y >= Ly) {
    pdata.y -= Ly;
  } else if (pdata.y < 0) {
    pdata.y += Ly;
  }
  moved = true;
}
