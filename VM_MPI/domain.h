/****************************************************************************/
/*           Decompose the simulation domain into several subdomains,       */
/*                       then run parallelly with MPI.                      */
/****************************************************************************/
#ifndef DOMAIN_H
#define DOMAIN_H
#include <vector>
#include <stack>
#include "mpi.h"
#include "rand.h"
#include "cell.h"

class BasicDomain
{
public:
  BasicDomain(double eta, double eps, double rho0,
              double Lx0, double Ly0, unsigned long long seed);
  ~BasicDomain();
  void create_particle_random(int nPar, double multiple);
  void create_from_snap(const std::string &filename);
  void output(int t);
  void create_cell_list();
  void update_position_edge_row(int row, double eta);
  void update_position_inner_rows(double eta);
  void pack(int row, double *buf, int &buf_size);
  void unpack(int row, const double *buf, int buf_size);
  void update_position(double eta);
  void send(int row, double *buf, int &buf_size, int dest, int tag,
            MPI_Request *req);
  void recv(double *buf, int buf_size, int source, int tag,
            MPI_Request *req);
  void accept(int row, double *buf, int *buf_size,
              MPI_Request *req, MPI_Status *stat);

protected:
  int tot_rank;
  int myrank;
  int pre_rank;
  int next_rank;

  double Lx;
  double Ly;
  double yl;
  double yh;

  Cell *cell;
  int ncols;
  int nrows;

  Node *particle;
  int MAX_PAR;
  int end_pos;
  int nPar_per_row;
  int nPar_per_task;
  std::stack <Node *> empty_pos;

  int MAX_BUF_SIZE;
  Ran *myran;

  std::ofstream fout_phi;
};

inline void BasicDomain::send(int row, double *buf, int &buf_size, int dest,
                              int tag, MPI_Request *req) {
  pack(row, buf, buf_size);
  MPI_Isend(buf, buf_size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, req);
}

inline void BasicDomain::recv(double *buf, int buf_size, int source,
                              int tag, MPI_Request *req) {
  MPI_Irecv(buf, buf_size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, req);
}

inline void BasicDomain::accept(int row, double *buf, int *buf_size,
                                MPI_Request *req, MPI_Status *stat) {
  MPI_Wait(req, stat);
  MPI_Get_count(stat, MPI_DOUBLE, buf_size);
  unpack(row, buf, *buf_size);
}

class StaticDomain :public BasicDomain
{
public:
  StaticDomain(double eta, double eps, double rho0,
               double Lx0, double Ly0, unsigned long long seed);
  ~StaticDomain();
  void update_velocity();
  void one_step(double eta, int t);
};

class DynamicDomain : public BasicDomain
{
public:
  DynamicDomain(double eta, double eps, double rho0,
                double Lx0, double Ly0, unsigned long long seed);
  ~DynamicDomain();
  void rearrange_domain(int t);
  void update_velocity();
  void one_step(double eta, int t);
private:
  Cell *cell_buf;
  int CELL_BUF_SIZE;
  int row_offset[2];
};

#endif
