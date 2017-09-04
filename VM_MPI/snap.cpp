#include "snap.h"

using namespace std;

Snap::Snap(const cmdline::parser & cmd, int tot_n):
           global_nPar(tot_n), idx_frame(0), offset0(0) {
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_rank);

  int nstep = cmd.get<int>("nstep");
  int dt = cmd.get<int>("cg_dt");
  double exp = cmd.get<double>("cg_exp");
  if (exp > 0 && dt > 0) {
    gene_log_frames(exp, cmd.get<int>("t_equil"));
    gene_lin_frames(dt, nstep, cmd.get<int>("t_equil"));
  } else if (exp > 0) {
    gene_log_frames(exp, nstep);
  } else if (dt > 0) {
    gene_lin_frames(dt, nstep, cmd.get<int>("t_equil"));
  }
  cout << "total frames: " << frame.size() << endl;
  dynamic = cmd.exist("dynamic") ? true : false;
}

Snap::~Snap() {
  MPI_File_close(&fh);
}

void Snap::gene_log_frames(double exponent, int t_end) {
  frame.push_back(1);
  double t = 1;
  int cur_frame = 1;
  while (true) {
    t *= exponent;
    if (t >= t_end) break;
    int tmp_frame = round(t);
    if (tmp_frame > cur_frame) {
      cur_frame = tmp_frame;
      frame.push_back(cur_frame);
    }
  }
}

void Snap::gene_lin_frames(int dt, int t_end, int t_equil) {
  for (int t = t_equil; t <= t_end; t++) {
    if (t % dt == 0)
      frame.push_back(t);
  }
}

void Snap::output(const Node * bird, int end_pos, int step) {
  if (step == frame[idx_frame]) {
    /* output some data */
    idx_frame++;
  }
}

void Snap::output_t_and_v(int t, double * sum_v) {
  double global_sum_v[2];
  MPI_Reduce(sum_v, global_sum_v, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myrank == 0) {
    global_sum_v[0] /= global_nPar;
    global_sum_v[1] /= global_nPar;
    MPI_File_write_at(fh, offset0, &t, 1, MPI_INT, MPI_STATUSES_IGNORE);
    MPI_File_write_at(fh, offset0 + sizeof(int), global_sum_v, 2,
                      MPI_DOUBLE, MPI_STATUSES_IGNORE);
  }
  offset0 += sizeof(int) + 2 * sizeof(double);
}

void Snap::set_offset(MPI_Offset *my_offset, MPI_Offset my_size) {
  if (!dynamic) {
    *my_offset = myrank * my_size + offset0;
  } else {
    MPI_Offset *size_v = new MPI_Offset[tot_rank];
    MPI_Offset *offset = new MPI_Offset[tot_rank];
    MPI_Gather(&my_size, 1, MPI_OFFSET,
               size_v, 1, MPI_OFFSET, 2, MPI_COMM_WORLD);
    if (myrank == 2) {
      MPI_Offset sum_size = 0;
      for (int i = 0; i < tot_rank; i++) {
        offset[i] = sum_size + offset0;
        sum_size += size_v[i];
      }
    }
    MPI_Scatter(offset, 1, MPI_OFFSET,
                my_offset, 1, MPI_OFFSET, 2, MPI_COMM_WORLD);
    delete[] size_v;
    delete[] offset;
  }
  offset0 += frame_size - 20;
}

CoarseGrainSnap::CoarseGrainSnap(const cmdline::parser &cmd, int tot_n):
                                 Snap(cmd, tot_n) {
  lx = cmd.get<int>("cg_lx");
  ly = cmd.get<int>("cg_ly");
  double Lx_global = cmd.get<double>("Lx");
  double Ly_global = cmd.get<double>("Ly");
  ncols = int(Lx_global / lx);
  nrows = int(Ly_global / ly / tot_rank);
  format = cmd.get<string>("cg_format");
  char fname[100];
  snprintf(fname, 100, "coarse%sc%s_%g_%g_%g_%g_%d_%d_%d_%llu.bin",
           delimiter.c_str(), format.c_str(),
           cmd.get<double>("eta"), cmd.get<double>("eps"),
           Lx_global, Ly_global, ncols, nrows * tot_rank, tot_n,
           cmd.get<unsigned long long>("seed"));
  cout << "file name " << fname << endl;
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                MPI_INFO_NULL, &fh);
  cout << "open " << fname << endl;
  ncells = ncols * nrows;
  global_ncells = ncells * tot_rank;
  if (format == "B") {
    frame_size = global_ncells * sizeof(unsigned char) + 20;
  } else if (format == "Hff") {
    frame_size = global_ncells * (sizeof(unsigned short) + 2 * sizeof(float)) + 20;
  }
  file_size = frame.size() * frame_size;
  MPI_File_preallocate(fh, file_size);
  cout << " size of file " << file_size << endl;
  cout << " size of one frame " << frame_size << endl;
}

CoarseGrainSnap::~CoarseGrainSnap() {
}

void CoarseGrainSnap::coarse_grain(const Node * bird, int end_pos, double yl,
                                   unsigned short * par_num,
                                   float * vx, float * vy, double *sum_v) {
  double y0 = yl + 1;
  sum_v[0] = sum_v[1] = 0;
  for (int i = 0; i < ncells; i++) {
    par_num[i] = 0;
    vx[i] = 0;
    vy[i] = 0;
  }
  for (int i = 0; i < end_pos; i++) {
    if (!bird[i].is_empty && !bird[i].is_ghost) {
      int col = bird[i].x / lx;
      if (col >= ncols) col = 0;
      int row = (bird[i].y - y0) / ly;
      if (row >= nrows) row = 0;
      int idx_cell = col + row * ncols;
      par_num[idx_cell]++;
      vx[idx_cell] += bird[i].vx;
      vy[idx_cell] += bird[i].vy;
      sum_v[0] += bird[i].vx;
      sum_v[1] += bird[i].vy;
    }
  }
  for (int i = 0; i < ncells; i++) {
    if (par_num[i] > 0) {
      vx[i] /= par_num[i];
      vy[i] /= par_num[i];
    }
  }
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
                            double yl, int nrows0, int t) {
  if (dynamic) {
    nrows = nrows0 - 2;
    ncells = nrows * ncols;
  }
  if (format == "B") {
    unsigned char *count = new unsigned char[ncells];
    double sum_v[2];
    coarse_grain_filtered(bird, end_pos, yl, count, sum_v[0], sum_v[1]);
    output_t_and_v(t, sum_v);
    MPI_Offset my_offset;
    set_offset(&my_offset, ncells * sizeof(unsigned char));
    MPI_File_write_at_all_begin(fh, my_offset, count, ncells, MPI_UNSIGNED_CHAR);
    MPI_File_write_at_all_end(fh, count, MPI_STATUS_IGNORE);
    delete[] count;
  } else if (format == "Hff") {
    unsigned short *particle_num = new unsigned short[ncells];
    float *vx = new float[ncells];
    float *vy = new float[ncells];
    double sum_v[2];
    coarse_grain(bird, end_pos, yl, particle_num, vx, vy, sum_v);
    // output time and vx, vy
    output_t_and_v(t, sum_v);
    // output particle number of each box
    MPI_Offset my_offset = offset0 + myrank * ncells * sizeof(unsigned short);
    MPI_File_write_at(fh, my_offset, particle_num, ncells,
                      MPI_UNSIGNED_SHORT, MPI_STATUSES_IGNORE);
    offset0 += tot_rank * ncells * sizeof(unsigned short);
    // output vx 
    my_offset = offset0 + myrank * ncells * sizeof(float);
    MPI_File_write_at(fh, my_offset, vx, ncells,
                      MPI_FLOAT, MPI_STATUSES_IGNORE);
    offset0 += tot_rank * ncells * sizeof(float);
    // output vy
    my_offset = offset0 + myrank * ncells * sizeof(float);
    MPI_File_write_at(fh, my_offset, vy, ncells,
                      MPI_FLOAT, MPI_STATUSES_IGNORE);
    offset0 += tot_rank * ncells * sizeof(float);
    delete[] particle_num;
    delete[] vx;
    delete[] vy;
  }
}

void CoarseGrainSnap::output(const Node * bird, int end_pos,
                             double yl, int nrows0, int t) {
  if (t == frame[idx_frame]) {
    write(bird, end_pos, yl, nrows0, t);
    if (myrank == 1)
      cout << idx_frame << "frame, t = " << t << endl;
    idx_frame++;
  }
}
