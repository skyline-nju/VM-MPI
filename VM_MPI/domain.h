#ifndef DOMAIN_H
#define DOMAIN_H
#include <vector>
#include <stack>
#include "mpi.h"
#include "rand.h"
#include "cell.h"

class SubDomain
{
public:
  SubDomain(double Lx_domain, double Ly_domain, int nrows_subdomain,
            int rank, unsigned long long seed);
  void create_particle_random(int nPar);
  void update_velocity_by_row(int row);
  void update_velocity_inner_rows();
  void update_position_inner_rows(double eta);
  void update_position_edge_row(double eta, int row);
  void create_cell_list();
  void remove_ghost_particle(int row);
  int count_valid_particle();
  void pack(int row, double *buff, int &buff_size);
  void unpack(int row, const double *buff, int buff_size);
  void update_velocity_MPI();
  void update_position_MPI(double eta);
  void one_step_MPI(double eta, int t);
private:
  double Lx;
  double Ly;
  double yl;
  double yh;
  int ncols;
  int nrows;
  int tot_rank;
  int myrank;
  int pre_rank;
  int next_rank;
  int MAX_BUFF_SIZE;
  Ran *myran;
  Cell *cell;
  std::vector <Node> particle;
  std::stack <unsigned int> empty_pos;
};
#endif
