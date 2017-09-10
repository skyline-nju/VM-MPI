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

  void set_first_row(double yl, int nrows_domain);

  template <typename T>
  void coarse_grain(T *par_num, double *sum_v, bool filtered,
                    const Node *bird, int end_pos) const;

  template <typename T1, typename T2>
  void coarse_grain(T1 *par_num, T2 *vx, T2 *vy, double *sum_v,
                    const Node *bird, int end_pos) const;
  template <typename T>
  void save_snap(T *par_num, MPI_Datatype type);
  
  template <typename T1, typename T2>
  void save_snap(T1 *par_num, MPI_Datatype type1,
                 T2 *vx, T2 *vy, MPI_Datatype type2);

  void output(const Node *bird, int end_pos, double yl,
              int nrows_domain, int t);

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

template<typename T>
void CoarseGrainSnap::coarse_grain(T * par_num, double * sum_v, bool filtered,
                                   const Node *bird, int end_pos) const {
  int ncells = ncols * nrows;
  for (int i = 0; i < ncells; i++) {
    par_num[i] = 0;
  }
  sum_v[0] = 0;
  sum_v[1] = 0;
  for (int i = 0; i < end_pos; i++) {
    if (!bird[i].is_empty && !bird[i].is_ghost) {
      sum_v[0] += bird[i].vx;
      sum_v[1] += bird[i].vy;
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

template<typename T1, typename T2>
void CoarseGrainSnap::coarse_grain(T1 *par_num, T2 *vx, T2 *vy, double *sum_v,
                                   const Node * bird, int end_pos) const {
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
      int row = bird[i].y / ly - first_row;
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

template<typename T>
void CoarseGrainSnap::save_snap(T * par_num, MPI_Datatype type) {
  MPI_Request req[2];
  T *recv_buf;
  if (flag_send) {
    MPI_Isend(par_num + (nrows - 1) * ncols, ncols, type,
              next_rank, 1, MPI_COMM_WORLD, &req[0]);
  }
  if (flag_recv) {
    recv_buf = new T[ncols];
    MPI_Irecv(recv_buf, ncols, type, pre_rank, 1, MPI_COMM_WORLD, &req[1]);
    MPI_Wait(&req[1], MPI_STATUS_IGNORE);
    for (int i = 0; i < ncols; i++) {
      par_num[i] += recv_buf[i];
    }
    delete[] recv_buf;
    MPI_Offset my_offset = offset0 + first_row * ncols * sizeof(T);
    MPI_File_write_at(fh, my_offset, par_num, (nrows - 1) * ncols,
                      type, MPI_STATUS_IGNORE);
  } else {
    MPI_Offset my_offset = offset0 + (first_row + 1) *  ncols * sizeof(T);
    MPI_File_write_at(fh, my_offset, par_num + ncols, (nrows - 2) * ncols,
                      type, MPI_STATUS_IGNORE);
  }
  offset0 += global_ncells * sizeof(T);
  if (flag_send)
    MPI_Wait(&req[0], MPI_STATUS_IGNORE);
}

template<typename T1, typename T2>
void CoarseGrainSnap::save_snap(T1 * par_num, MPI_Datatype type1,
                                T2 * vx, T2 * vy, MPI_Datatype type2) {
  MPI_Request req[6];
  T1 *recv_buf0;
  T2 *recv_buf1;
  T2 *recv_buf2;
  if (flag_send) {
    int k = (nrows - 1) * ncols;
    MPI_Isend(&par_num[k], ncols, type1, next_rank, 0, MPI_COMM_WORLD, &req[0]);
    MPI_Isend(&vx[k], ncols, type2, next_rank, 1, MPI_COMM_WORLD, &req[1]);
    MPI_Isend(&vy[k], ncols, type2, next_rank, 2, MPI_COMM_WORLD, &req[2]);
  }
  if (flag_recv) {
    recv_buf0 = new T1[ncols];
    recv_buf1 = new T2[ncols];
    recv_buf2 = new T2[ncols];
    MPI_Irecv(recv_buf0, ncols, type1, pre_rank, 0, MPI_COMM_WORLD, &req[3]);
    MPI_Irecv(recv_buf1, ncols, type2, pre_rank, 1, MPI_COMM_WORLD, &req[4]);
    MPI_Irecv(recv_buf2, ncols, type2, pre_rank, 2, MPI_COMM_WORLD, &req[5]);
    for (int i = ncols; i < (nrows - 1) * ncols; i++) {
      if (par_num[i] != 0) {
        vx[i] /= par_num[i];
        vy[i] /= par_num[i];
      }
    }
    MPI_Waitall(3, &req[3], MPI_STATUSES_IGNORE);
    for (int i = 0; i < ncols; i++) {
      par_num[i] += recv_buf0[i];
      vx[i] += recv_buf1[i];
      vy[i] += recv_buf2[i];
      if (par_num[i] != 0) {
        vx[i] /= par_num[i];
        vy[i] /= par_num[i];
      }
    }
    delete[] recv_buf0;
    delete[] recv_buf1;
    delete[] recv_buf2;

    MPI_Offset my_offset1 = offset0 + first_row * ncols * sizeof(T1);
    offset0 += global_ncells * sizeof(T1);
    MPI_Offset my_offset2 = offset0 + first_row * ncols * sizeof(T2);
    offset0 += global_ncells * sizeof(T2);
    MPI_Offset my_offset3 = offset0 + first_row * ncols * sizeof(T2);
    offset0 += global_ncells * sizeof(T2);
    int count = (nrows - 1) * ncols;

    MPI_File_write_at(fh, my_offset1, par_num, count, type1, MPI_STATUS_IGNORE);
    MPI_File_write_at(fh, my_offset2, vx, count, type2, MPI_STATUS_IGNORE);
    MPI_File_write_at(fh, my_offset3, vy, count, type2, MPI_STATUS_IGNORE);
  } else {
    MPI_Offset my_offset1 = offset0 + (first_row + 1) * ncols * sizeof(T1);
    offset0 += global_ncells * sizeof(T1);
    MPI_Offset my_offset2 = offset0 + (first_row + 1) * ncols * sizeof(T2);
    offset0 += global_ncells * sizeof(T2);
    MPI_Offset my_offset3 = offset0 + (first_row + 1) * ncols * sizeof(T2);
    offset0 += global_ncells * sizeof(T2);
    int count = (nrows - 2) * ncols;

    MPI_Request req;
    MPI_File_iwrite_at(fh, my_offset1, &par_num[ncols], count, type1, &req);
    for (int i = ncols; i < (nrows - 1) * ncols; i++) {
      if (par_num[i] != 0) {
        vx[i] /= par_num[i];
        vy[i] /= par_num[i];
      }
    }
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    MPI_File_write_at(fh, my_offset2, &vx[ncols], count, type2, MPI_STATUS_IGNORE);
    MPI_File_write_at(fh, my_offset3, &vy[ncols], count, type2, MPI_STATUS_IGNORE);
  }
  if (flag_send) {
    MPI_Waitall(3, &req[0], MPI_STATUSES_IGNORE);
  }
}

#endif

