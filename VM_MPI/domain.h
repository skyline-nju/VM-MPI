#ifndef DOMAIN_H
#define DOMAIN_H
#include <vector>
#include <stack>
#include "mpi.h"
#include "rand.h"
#include "cell.h"

class StaticDomain
{
public:
  StaticDomain(double Lx_domain, double Ly_domain, int ntask, int rank,
               unsigned long long seed, double eta, double eps, double rho0);
  ~StaticDomain();
  void create_particle_random(int nPar);
  void create_from_snap(const std::string filename);
  void comm_start(int src_row, int &dest_row, double **sbuff, double **rbuff,
                  MPI_Request *sreq, MPI_Request *rreq);
  void comm_end(int dest_row, double *sbuff, double *rbuff,
                MPI_Request *sreq, MPI_Request *rreq);
  void update_velocity_MPI();
  void update_position_MPI(double eta);
  void one_step_MPI(double eta, int t);
  void output(int t);
private:
  double Lx;
  double Ly;
  double yl;
  double yh;

  int tot_rank;
  int myrank;
  int pre_rank;
  int next_rank;

  int ncols;
  int nrows;
  Cell *cell;

  Node *particle;
  int MAX_PAR;
  int end_pos;
  std::stack <Node *> empty_pos;

  int MAX_BUFF_SIZE;
  Ran *myran;

  std::ofstream fout_phi;
};

class DynamicDomain
{
public:
  DynamicDomain(double Lx_domain, double Ly_domain, int ntask, int rank,
                unsigned long long seed, double eta, double eps, double rho0);
  ~DynamicDomain();
  void create_particle_random(int nPar);
  void create_from_snap(const std::string filename);
  void rearrange();
  //void comm_start(int src_row, int &dest_row, double **sbuff, double **rbuff,
  //  MPI_Request *sreq, MPI_Request *rreq);
  //void comm_end(int dest_row, double *sbuff, double *rbuff,
  //  MPI_Request *sreq, MPI_Request *rreq);
  void update_velocity();
  void update_position(double eta);
  //void one_step_MPI(double eta, int t);
  //void output(int t);
private:
  int tot_rank;
  int myrank;
  int pre_rank;
  int next_rank;

  double Lx;
  double Ly;
  double yl;
  double yh;

  Cell *cell_buff;
  Cell *cell;
  int MAX_CELL;      /* Max number of cell */
  int ncols;
  int nrows;
  int row_offset[2];

  Node *particle;
  int MAX_PAR_BUFF;       /* Max size for particle buff */
  int end_pos;       /* Last index */
  int nPar_per_row;
  int nPar_per_task;
  std::stack <Node *> empty_pos;

  int MAX_BUFF_SIZE;

  Ran *myran;
  std::ofstream fout_phi;
};
#endif
