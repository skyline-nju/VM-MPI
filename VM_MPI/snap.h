#ifndef SNAP_H
#define SNAP_H
#include "comn.h"
#include "node.h"
#include "cell.h"
#include "mpi.h"
#include "cmdline.h"

class Snap
{
public:
  Snap(const cmdline::parser &cmd, int tot_n);
  ~Snap();
  void gene_log_frames(double exponent, int t_end);
  void gene_lin_frames(int dt, int t_end, int t_equil);
  void output(const Node *bird, int end_pos, int step);
  void output_t_and_v(int t, double *sum_v);
  void set_offset(MPI_Offset *my_offset, MPI_Offset my_size);

protected:
  MPI_File fh;
  std::vector<int> frame;
  MPI_Offset frame_size;
  MPI_Offset file_size;
  MPI_Offset offset0;
  int idx_frame;
  int global_nPar;
  int tot_rank;
  int myrank;
  int pre_rank;
  int next_rank;
  bool dynamic;
};

class CoarseGrainSnap : public Snap
{
public:
  CoarseGrainSnap(const cmdline::parser &cmd, int tot_n);
  ~CoarseGrainSnap();
  void output(const Node *bird, int end_pos, double yl,
              int nrows_domain, int t);
  void set_first_row(double yl, int nrows_domain);
  void save_snap(unsigned char *par_num);
  void save_snap(unsigned short *par_num, float * vx, float * vy);

private:
  bool flag_recv;
  bool flag_send;
  int lx;
  int ly;
  int ncols;
  int nrows;
  int ncells;
  int first_row;
  int global_ncells;
  std::string format;
};

template <typename T>
void coarse_grain(const Node *bird, int end_pos,
                  int first_row, int lx, int ly, int ncols, int nrows,
                  T *par_num, double *sum_v, bool filtered) {
  int ncells = ncols * nrows;
  for (int i = 0; i < ncells; i++) {
    par_num[i] = 0;
  }
  sum_v[0] = 0;
  sum_v[1] = 0;
  for (int i = 0; i < end_pos; i++) {
    if (!bird[i].is_empty && !bird[i].is_ghost) {
      svx += bird[i].vx;
      svy += bird[i].vy;
      if (!filtered || bird[i].vx > 0) {
        int col = bird[i].x / lx;
        if (col > ncols) col = 0;
        int row = bird[i].y / ly - first_row;
        if (row > nrows) row = 0;
        int idx_cell = col + row * ncols;
        par_num[idx_cell]++;
      }
    }
  }
}

template <typename T1, typename T2>
void coarse_grain(const Node *bird, int end_pos,
                  int first_row, int lx, int ly,
                  int ncols, int nrows,
                  T1 * par_num, T2 * vx, T2 * vy, double * sum_v) {
  int ncells = ncols * nrows;
  for (int i = 0; i < ncells; i++) {
    par_num[i] = 0;
    vx[i] = vy[i] = 0;
  }
  sum_v[0] = 0;
  sum_v[1] = 0;
  for (int i = 0; i < end_pos; i++) {
    if (!bird[i].is_empty && !bird[i].is_ghost) {
      int col = bird[i].x / lx;
      if (col >= ncols) col = 0;
      int row = bird[i].y  / ly - first_row;
      if (row >= nrows) row = 0;
      int idx_cell = col + row * ncols;
      par_num[idx_cell]++;
      vx[idx_cell] += bird[i].vx;
      vy[idx_cell] += bird[i].vy;
      sum_v[0] += bird[i].vx;
      sum_v[1] += bird[i].vy;
    }
  }
}

#endif
