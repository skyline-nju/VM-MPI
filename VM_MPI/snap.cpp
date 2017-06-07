#include "snap.h"

using namespace std;

Snap::Snap(const cmdline::parser & cmd, int tot_n):
           global_nPar(tot_n), idx(0), offset0(0) {
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_rank);
  char str[100];
  snprintf(str, 100, "c_%g_%g_%g_%g_%d_%llu",
           cmd.get<double>("eta"), cmd.get<double>("eps"),
           cmd.get<double>("Lx"), cmd.get<double>("Ly"), tot_n,
           cmd.get<unsigned long long>("seed"));
  int nstep = cmd.get<int>("nstep");
  int dt = cmd.get<int>("cg_dt");
  double exp = cmd.get<double>("cg_exp");
  dynamic = cmd.exist("dynamic") ? true : false;
  char f1[100];
  char f2[100];
  if (dt > 0) {
    snprintf(f1, 100, "%s_%d.bin", str, dt);
    snprintf(f2, 100, "%s_%d.dat", str, dt);
    gene_lin_frames(dt, nstep);
  } else {
    snprintf(f1, 100, "%s_%g.bin", str, exp);
    snprintf(f2, 100, "%s_%g.dat", str, exp);
    gene_log_frames(exp, nstep);
  }
  MPI_File_open(MPI_COMM_WORLD, f1, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                MPI_INFO_NULL, &fh);
  if (myrank == 2) {
    fout.open(f2);
  }
}

Snap::~Snap() {
  MPI_File_close(&fh);
  if (myrank == 2)
    fout.close();
}

void Snap::gene_log_frames(double exponent, int t_end) {
  frame.push_back(1);
  double t = 1;
  int cur_frame = 1;
  while (true) {
    t *= exponent;
    if (t > t_end) break;
    int tmp_frame = round(t);
    if (tmp_frame > cur_frame) {
      cur_frame = tmp_frame;
      frame.push_back(cur_frame);
    }
  }
}

void Snap::gene_lin_frames(int dt, int t_end) {
  for (int t = 1; t <= t_end; t++) {
    if (t % dt == 0)
      frame.push_back(t);
  }
}

void Snap::output(const Node * bird, int end_pos, int step) {
  if (step == frame[idx]) {
    /* output some data */
    idx++;
  }
}

void Snap::output_phi(const double *sv) {
  double g_sum_v[2];
  MPI_Reduce(sv, g_sum_v, 2, MPI_DOUBLE, MPI_SUM, 2, MPI_COMM_WORLD);
  if (myrank == 2) {
    double phi = sqrt(g_sum_v[0] * g_sum_v[0] +
                 g_sum_v[1] * g_sum_v[1]) / global_nPar;
    fout << frame[idx] << "\t" << phi << endl;
  }
}

void Snap::set_offset(MPI_Offset *my_offset, MPI_Offset my_size) {
  if (!dynamic) {
    *my_offset = myrank * my_size + offset0;
  } else {
    int *size_v = new int[tot_rank];
    MPI_Offset *offset = new MPI_Offset[tot_rank];
    MPI_Gather(&my_size, 1, MPI_OFFSET,
               size_v, 1, MPI_OFFSET, 2, MPI_COMM_WORLD);
    if (myrank == 2) {
      MPI_Offset sum_size = 0;
      for (int i = 0; i < tot_rank; i++) {
        offset[i] = sum_size + offset0;
        sum_size += size_v[i];
        cout << "offset = " << offset0 << endl;
      }
    }
    MPI_Scatter(offset, 1, MPI_OFFSET,
                my_offset, 1, MPI_OFFSET, 2, MPI_COMM_WORLD);
    delete[] size_v;
    delete[] offset;
  }
  offset0 += frame_size;
}

CoarseGrainSnap::CoarseGrainSnap(const cmdline::parser &cmd, int tot_n):
                                 Snap(cmd, tot_n) {
  ncols = int(cmd.get<double>("Lx"));
  nrows = int(cmd.get<double>("Ly") / tot_rank);
  ncells = ncols * nrows;
  global_ncells = ncells * tot_rank;
  frame_size = global_ncells * sizeof(unsigned char);
  file_size = frame.size() * frame_size;
  MPI_File_preallocate(fh, file_size);
  cout << " size of file " << file_size << endl;
}

CoarseGrainSnap::~CoarseGrainSnap() {
}

void CoarseGrainSnap::coarse_grain_filtered(const Node * bird, int end_pos,
                                            double yl, unsigned char *count,
                                            double &svx, double &svy) {
  double y0 = yl + 1;
  svx = svy = 0;
  for (int i = 0; i < ncells; i++)
    count[i] = 0;
  for (int i = 0; i < end_pos; i++) {
    if (!bird[i].is_empty && !bird[i].is_ghost) {
      svx += bird[i].vx;
      svy += bird[i].vy;
      if (bird[i].vx > 0) {
        int col = int(bird[i].x);
        int row = int(bird[i].y - y0);
        if (col + ncols * row >= ncells) {
          cout << "col = " << col << "\trow = " << row << "\tncells = " << ncells << endl;
        }
        count[col + ncols * row]++;
      }
    }
  }
}

void CoarseGrainSnap::write(const Node *bird, int end_pos,
                            double yl, int nrows0) {
  nrows = nrows0 - 2;
  ncells = nrows * ncols;
  unsigned char *count = new unsigned char[ncells];
  double sum_v[2];
  coarse_grain_filtered(bird, end_pos, yl, count, sum_v[0], sum_v[1]);
  MPI_Offset my_offset;
  set_offset(&my_offset, ncells * sizeof(unsigned char));
  MPI_File_write_at_all_begin(fh, my_offset, count, ncells, MPI_UNSIGNED_CHAR);
  output_phi(sum_v);
  MPI_File_write_at_all_end(fh, count, MPI_STATUS_IGNORE);
  delete[] count;
}

void CoarseGrainSnap::output(const Node * bird, int end_pos,
                             double yl, int nrows0, int t) {
  if (t == frame[idx]) {
    write(bird, end_pos, yl, nrows0);
    idx++;
  }
}
