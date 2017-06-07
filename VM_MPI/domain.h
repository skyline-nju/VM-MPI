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
#include "cmdline.h"
#include "snap.h"

class BasicDomain
{
public:
  BasicDomain(double eta, double eps, double rho0,
              double Lx0, double Ly0, unsigned long long seed);
  BasicDomain(const cmdline::parser &cmd);
  virtual ~BasicDomain() {};
  void create_particle_random(int tot_nPar, double magnification);
  void create_from_snap(const std::string &filename, double magnification);
  void create_from_snap(const cmdline::parser &cmd, double magnification);
  void output(int t);
  void ini_output(const std::string &tag, double eta, double eps,
                  double rho0, unsigned long long seed);
  void create_cell_list();
  void ini_disorder(double eps, double *disorder, int n);
  void update_position_edge_row(int row, double eta);
  void update_position_inner_rows(double eta);
  void pack(int row, double *buf, int &buf_size);
  void unpack(int row, const double *buf, int buf_size);
  void send(int row, double *buf, int &buf_size, int dest, int tag,
            MPI_Request *req);
  void recv(double *buf, int buf_size, int source, int tag,
            MPI_Request *req);
  void accept(int row, double *buf, int *buf_size,
              MPI_Request *req, MPI_Status *stat);
  double get_noise(double eta);
  void update_position(double eta);
  void update_velocity();
  void shift_boundary(int pre_nPar, int next_nPar, int *offset);
  virtual void one_step(double eta, int t) = 0;

protected:
  int tot_rank;
  int myrank;
  int pre_rank;
  int next_rank;

  double eta;
  double eps;
  double rho0;
  double Lx;
  double Ly;
  double yl;
  double yh;
  int tot_steps;

  Cell *cell;
  int ncols;
  int nrows;

  Node *particle;
  int my_nPar;
  int MAX_PAR;
  int end_pos;
  int nPar_per_row;
  int nPar_per_task;
  std::stack <Node *> empty_pos;

  int MAX_BUF_SIZE;
  Ran *myran;
  bool disorder_free;
  bool flag_coarse_grain;

  std::time_t beg_time;

  CoarseGrainSnap *cgs;       // output coarse-grained snapshot
  std::ofstream fout;         // output order parameters for rank = 0
                              // output log for rank = 1

};


class StaticDomain :public BasicDomain
{
public:
  StaticDomain(double eta, double eps, double rho0,
               double Lx0, double Ly0, unsigned long long seed);
  StaticDomain(const cmdline::parser &cmd);
  ~StaticDomain();
  void set_cell_list(double eps);
  void one_step(double eta, int t);
};

class DynamicDomain : public BasicDomain
{
public:
  DynamicDomain(double eta, double eps, double rho0,
                double Lx0, double Ly0, unsigned long long seed,
                int refresh_rate0);
  DynamicDomain(const cmdline::parser &cmd);
  ~DynamicDomain();

  void set_cell_list(double eps);
  void global_rearrange(int t);       /* update row_offset */
  void local_rearrange(int t);
  void update_velocity_dynamic();
  void one_step(double eta, int t);

private:
  Cell *cell_buf;     // The size of cell_buf is larger than the size of
  int CELL_BUF_SIZE;  // cell, so that there is enough space for cell to
  int MAX_ROW;        // expand or shrink.
                      
  int refresh_rate;
  int row_offset[2];  // Control the size of subdomain. row_offset[0] = 1
                      // means the bottom boundary of the subdomain moves
                      // downward by 1, while row_offset[1] = 1 means the
                      // top boundary moves upward by 1.
};

inline double BasicDomain::get_noise(double eta) {
  return eta * 2.0 * PI * (myran->doub() - 0.5);
}
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
#endif
