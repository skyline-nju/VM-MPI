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
  bool dynamic;
};

class CoarseGrainSnap : public Snap
{
public:
  CoarseGrainSnap(const cmdline::parser &cmd, int tot_n);
  ~CoarseGrainSnap();
  void coarse_grain(const Node *bird, int end_pos, double yl,
                    unsigned short *par_num, float *vx, float *vy,
                    double *sum_v);
  void coarse_grain_filtered(const Node *bird, int end_pos, double yl, 
                             unsigned char *count, double &svx, double &svy);
  void write(const Node *bird, int end_pos, double yl, int nrows0, int t);
  void output(const Node *bird, int end_pos, double yl, int nrows0, int t);
private:
  int lx;
  int ly;
  int ncols;
  int nrows;
  int ncells;
  int global_ncells;
  std::string format;
};
#endif
