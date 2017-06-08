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
  void gene_lin_frames(int dt, int t_end);
  void output(const Node *bird, int end_pos, int step);
  void output_phi(const double *sv);
  void set_offset(MPI_Offset *my_offset, MPI_Offset my_size);

protected:
  MPI_File fh;
  std::ofstream fout;
  std::vector<int> frame;
  MPI_Offset frame_size;
  MPI_Offset file_size;
  MPI_Offset offset0;
  int idx;
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
  void coarse_grain_filtered(const Node *bird, int end_pos, double yl, 
                             unsigned char *count, double &svx, double &svy);
  void write(const Node *bird, int end_pos, double yl, int nrows0);
  void output(const Node *bird, int end_pos, double yl, int nrows0, int t);
private:
  int ncols;
  int nrows;
  int ncells;
  int global_ncells;
};
#endif
