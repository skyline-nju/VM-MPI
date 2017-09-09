#include "snap.h"

using namespace std;

Snap::Snap(const cmdline::parser & cmd, int tot_n):
           global_nPar(tot_n), idx_frame(0), offset0(0) {
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_rank);
  pre_rank = myrank == 0 ? tot_rank - 1 : myrank - 1;
  next_rank = myrank == tot_rank - 1 ? 0 : myrank + 1;

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
  double gl_Lx = cmd.get<double>("Lx");
  double gl_Ly = cmd.get<double>("Ly");
  int gl_ncols = int(gl_Lx / lx);
  int gl_nrows = int(gl_Ly / ly);
  ncols = gl_ncols;
  format = cmd.get<string>("cg_format");
  char fname[100];
  snprintf(fname, 100, "coarse%sc%s_%g_%g_%g_%g_%d_%d_%d_%llu.bin",
           delimiter.c_str(), format.c_str(),
           cmd.get<double>("eta"), cmd.get<double>("eps"),
           gl_Lx, gl_Ly, gl_ncols, gl_nrows, tot_n,
           cmd.get<unsigned long long>("seed"));
  cout << "file name " << fname << endl;
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                MPI_INFO_NULL, &fh);
  cout << "open " << fname << endl;
  global_ncells = gl_ncols * gl_nrows;
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


void CoarseGrainSnap::set_first_row(double yl, int nrows_domain) {
  int ymin = yl + 1;
  int ymax = ymin + nrows_domain - 2;
  cout << "ymin = " << ymin << "\tymax = " << ymax << endl;
  if (ymin % ly == 0) {
    flag_recv = false;
    first_row = ymin / ly - 1;
    int dy = ymax - ymin + ly;
    flag_send = dy % ly == 0 ? false : true;
    nrows = dy / ly + 1;
  } else {
    flag_recv = true;
    first_row = ymin / ly;
    int dy = ymax - ymin + ymin % ly;
    flag_send = dy % ly == 0 ? false : true;
    nrows = dy / ly + 1;
  }
  ncells = nrows * ncols;
}

void CoarseGrainSnap::save_snap(unsigned char * par_num) {
  MPI_Request req[2];
  unsigned char *recv_buf;
  if (flag_send) {
    MPI_Isend(par_num + (nrows - 1) * ncols, ncols, MPI_UNSIGNED_CHAR,
              next_rank, 1, MPI_COMM_WORLD, &req[0]);
  }
  if (flag_recv) {
    recv_buf = new unsigned char[ncols];
    MPI_Irecv(recv_buf, ncols, MPI_UNSIGNED_CHAR,
              pre_rank, 1, MPI_COMM_WORLD, &req[1]);
    MPI_Wait(&req[1], MPI_STATUS_IGNORE);
    for (int i = 0; i < ncols; i++) {
      par_num[i] += recv_buf[i];
    }
    delete[] recv_buf;
    MPI_Offset my_offset = offset0 + first_row * ncols * sizeof(unsigned char);
    MPI_File_write_at(fh, my_offset, par_num, (nrows - 1) * ncols,
                      MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
  } else {
    MPI_Offset my_offset = offset0 +  (first_row + 1) *  ncols * sizeof(unsigned char);
    MPI_File_write_at(fh, my_offset, par_num + ncols, (nrows - 2) * ncols,
                      MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
  }
  offset0 += global_ncells * sizeof(unsigned char);
  if (flag_send)
    MPI_Wait(&req[0], MPI_STATUS_IGNORE);
}

void CoarseGrainSnap::save_snap(unsigned short *par_num, float *vx, float *vy) {
  MPI_Request req[6];
  unsigned short *recv_buf0;
  float *recv_buf1;
  float *recv_buf2;
  if (flag_send) {
    int k = (nrows - 1) * ncols;
    MPI_Isend(&par_num[k], ncols, MPI_UNSIGNED_SHORT, next_rank, 0,
              MPI_COMM_WORLD, &req[0]);
    MPI_Isend(&vx[k], ncols, MPI_FLOAT, next_rank, 1, MPI_COMM_WORLD, &req[1]);
    MPI_Isend(&vy[k], ncols, MPI_FLOAT, next_rank, 2, MPI_COMM_WORLD, &req[2]);
  }
  if (flag_recv) {
    recv_buf0 = new unsigned short[ncols];
    recv_buf1 = new float[ncols];
    recv_buf2 = new float[ncols];
    MPI_Irecv(recv_buf0, ncols, MPI_UNSIGNED_SHORT, pre_rank, 0,
              MPI_COMM_WORLD, &req[3]);
    MPI_Irecv(recv_buf1, ncols, MPI_FLOAT, pre_rank, 1,
              MPI_COMM_WORLD, &req[4]);
    MPI_Irecv(recv_buf2, ncols, MPI_FLOAT, pre_rank, 2,
              MPI_COMM_WORLD, &req[5]);
    for (int i = ncols; i < (nrows - 1) * ncols; i++) {
      if (par_num[i] != 0) {
        vx[i] /= par_num[i];
        vy[i] /= par_num[i];
      }
    }
    MPI_Waitall(3, &req[3], MPI_STATUSES_IGNORE);
    for (int i = 0; i < ncols; i++) {
      par_num[i] += recv_buf0[i];
      if (par_num[i] != 0) {
        vx[i] = (vx[i] + recv_buf1[i]) / par_num[i];
        vy[i] = (vy[i] + recv_buf2[i]) / par_num[i];
      }
    }
    int sum = 0;
    for (int i = 0; i < (nrows-1) * ncols; i++) {
      sum += par_num[i];
    }
    cout << myrank << "\tnum = " << sum << endl;
    delete[] recv_buf0;
    delete[] recv_buf1;
    delete[] recv_buf2;
    MPI_Offset my_offset = offset0 + first_row * ncols * sizeof(unsigned short);
    MPI_File_write_at(fh, my_offset, par_num, (nrows - 1) * ncols,
                      MPI_UNSIGNED_SHORT, MPI_STATUS_IGNORE);
    offset0 += global_ncells * sizeof(unsigned short);
    my_offset = offset0 + first_row * ncols * sizeof(float);
    MPI_File_write_at(fh, my_offset, vx, (nrows - 1) * ncols,
                      MPI_FLOAT, MPI_STATUS_IGNORE);
    offset0 += global_ncells * sizeof(float);
    my_offset = offset0 + first_row * ncols * sizeof(float);
    MPI_File_write_at(fh, my_offset, vy, (nrows - 1) * ncols,
                      MPI_FLOAT, MPI_STATUS_IGNORE);
    offset0 += global_ncells * sizeof(float);
  } else {
    MPI_Offset my_offset = offset0 + (first_row + 1) * ncols * sizeof(unsigned short);
    MPI_Request req;
    MPI_File_iwrite_at(fh, my_offset, &par_num[ncols], (nrows - 2) * ncols,
                       MPI_UNSIGNED_SHORT, &req);
    for (int i = ncols; i < (nrows - 1) * ncols; i++) {
      if (par_num[i] != 0) {
        vx[i] /= par_num[i];
        vy[i] /= par_num[i];
      }
    }
    int sum = 0;
    for (int i = ncols; i < (nrows - 1) * ncols; i++) {
      sum += par_num[i];
    }
    cout << myrank << "\tnum = " << 512 * 512 - sum << endl;
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    offset0 += global_ncells * sizeof(unsigned short);
    my_offset = offset0 + (first_row + 1) * ncols * sizeof(float);
    MPI_File_write_at(fh, my_offset, &vx[ncols], (nrows - 2) * ncols,
                      MPI_FLOAT, MPI_STATUS_IGNORE);
    offset0 += global_ncells * sizeof(float);
    my_offset = offset0 + (first_row + 1) * ncols * sizeof(float);
    MPI_File_write_at(fh, my_offset, &vy[ncols], (nrows - 2) * ncols,
                      MPI_FLOAT, MPI_STATUS_IGNORE);
    offset0 += global_ncells * sizeof(float);
  }
  if (flag_send) {
    MPI_Waitall(3, &req[0], MPI_STATUSES_IGNORE);
  }
}

void CoarseGrainSnap::output(const Node * bird, int end_pos,
                             double yl, int nrows_domain, int t) {
  if (t == frame[idx_frame]) {
    set_first_row(yl, nrows_domain);
    double sum_v[2];
    if (format == "Hff") {
      unsigned short *par_num = new unsigned short[ncells];
      float *vx = new float[ncells];
      float *vy = new float[ncells];
      coarse_grain(bird, end_pos, first_row, lx, ly, ncols, nrows,
                   par_num, vx, vy, sum_v);
      output_t_and_v(t, sum_v);
      save_snap(par_num, vx, vy);
      delete[] par_num;
      delete[] vx;
      delete[] vy;
    } else {
      unsigned char *par_num = new unsigned char[ncells];
      coarse_grain(bird, end_pos, first_row, lx, ly, ncols, nrows,
                   par_num, sum_v, false);
      int sum = 0;
      for (int i = ncols; i < (nrows - 1) * ncols; i++)
        sum += par_num[i];
      output_t_and_v(t, sum_v);
      save_snap(par_num);
      delete[] par_num;
    }
    if (myrank == 1)
      cout << idx_frame << "frame, t = " << t << endl;
    idx_frame++;
  }
}
