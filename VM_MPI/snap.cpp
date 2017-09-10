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
  MPI_File_set_size(fh, file_size);
  MPI_File_close(&fh);
}


void CoarseGrainSnap::set_first_row(double yl, int nrows_domain) {
  int ymin = yl + 1;
  int ymax = ymin + nrows_domain - 2;
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

void CoarseGrainSnap::output(const Node * bird, int end_pos,
                             double yl, int nrows_domain, int t) {
  if (t == frame[idx_frame]) {
    set_first_row(yl, nrows_domain);
    double sum_v[2];
    if (format == "Hff") {
      unsigned short *par_num = new unsigned short[ncells];
      float *vx = new float[ncells];
      float *vy = new float[ncells];
      coarse_grain(par_num, vx, vy, sum_v, bird, end_pos);
      output_t_and_v(t, sum_v);
      save_snap(par_num, MPI_UNSIGNED_SHORT, vx, vy, MPI_FLOAT);
      delete[] par_num;
      delete[] vx;
      delete[] vy;
    } else if (format == "B") {
      unsigned char *par_num = new unsigned char[ncells];
      coarse_grain(par_num, sum_v, false, bird, end_pos);
      output_t_and_v(t, sum_v);
      save_snap(par_num, MPI_UNSIGNED_CHAR);
      delete[] par_num;
    }
    cout << myrank << "\t" << t << "\t" << offset0 << endl;
    idx_frame++;
  }
}
