#ifndef DOMAIN_H
#define DOMAIN_H
#include "rand.h"
#include "cell.h"
#include <vector>
#include <stack>

class SubDomain
{
public:
  SubDomain(double _Lx, double _Ly, int tot_domain, int idx,
            unsigned long long seed);
  void update_velocity_inner();
  void update_velocity_edge();
  void update_position_inner(double eta);
  void update_position_edge(double eta);
  int get_pNum(int row);
  void pack(int row, std::vector<ParticleData> &buff);
  void unpack(const std::vector<ParticleData> &buff, int row);
private:
  double Lx;
  double Ly;
  double Ly_l;
  double Ly_h;
  int ncols;
  int nrows;
  Cell *cell;
  std::vector <Node> particle;
  std::stack <unsigned int> empty_pos;

  int tot_cells;
  Ran *myran;
};
#endif
