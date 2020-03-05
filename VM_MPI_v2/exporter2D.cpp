#include "exporter2D.h"

void exporter::ExporterBase::set_lin_frame(int start, int n_step, int sep) {
  n_step_ = n_step;
  for (auto i = start + sep; i <= n_step_; i += sep) {
    frames_arr_.push_back(i);
  }
  frame_iter_ = frames_arr_.begin();
}

bool exporter::ExporterBase::need_export(int i_step) {
  bool flag = false;
  if (!frames_arr_.empty() && i_step == (*frame_iter_)) {
    frame_iter_++;
    flag = true;
  }
  return flag;
}

exporter::LogExporter::LogExporter(const std::string& outfile, int start, int n_step, int sep, int np, MPI_Comm group_comm)
  : ExporterBase(start, n_step, sep), n_par_(np), comm_(group_comm) {
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(comm_, &my_rank);
  if (my_rank == 0) {
#endif
    fout.open(outfile);
    t_start_ = std::chrono::system_clock::now();
    auto start_time = std::chrono::system_clock::to_time_t(t_start_);
    char str[100];
    std::strftime(str, 100, "%c", std::localtime(&start_time));
    fout << "Started simulation at " << str << "\n";
#ifdef USE_MPI
  }
#endif
}

exporter::LogExporter::~LogExporter() {
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(comm_, &my_rank);
  if (my_rank == 0) {
#endif
    const auto t_now = std::chrono::system_clock::now();
    auto end_time = std::chrono::system_clock::to_time_t(t_now);
    char str[100];
    // ReSharper disable CppDeprecatedEntity
    std::strftime(str, 100, "%c", std::localtime(&end_time));
    // ReSharper restore CppDeprecatedEntity
    fout << "Finished simulation at " << str << "\n";
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    fout << "speed=" << std::scientific << step_count_ * double(n_par_) / elapsed_seconds.count()
      << " particle time step per second per core\n";
    fout.close();
#ifdef USE_MPI
  }
#endif
}

void exporter::LogExporter::record(int i_step) {
  bool flag;
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(comm_, &my_rank);
  flag = my_rank == 0 && need_export(i_step);
#else
  flag = need_export(i_step);
#endif
  if (flag) {
    const auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    const auto dt = elapsed_seconds.count();
    const auto hour = int(dt / 3600);
    const auto min = int((dt - hour * 3600) / 60);
    const int sec = dt - hour * 3600 - min * 60;
    fout << i_step << "\t" << hour << ":" << min << ":" << sec << std::endl;
  }
  step_count_++;
}

exporter::OrderParaExporter_2::OrderParaExporter_2(const std::string& outfile, int start, int n_step, int sep,
                                                   const Vec_2<double>& gl_l, MPI_Comm group_comm)
  : ExporterBase(start, n_step, sep), comm_(group_comm) {
  int tot_proc = 1;
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(comm_, &my_rank);
  MPI_Comm_size(comm_, &tot_proc);
  if (my_rank == 0) {
#endif
    fout_.open(outfile);
#ifdef USE_MPI
  }
#endif
  if (gl_l.x == gl_l.y && (int(gl_l.x) % 32) == 0) {
    flag_phi_box_ = true;
    for (int i = 32; i <= int(gl_l.x); i *= 2) {
      L_arr_.push_back(i);
    }
    max_cells_ = int(int(gl_l.x * gl_l.y) / (32 * 32) / tot_proc * 1.5);
  } else {
    flag_phi_box_ = false;
    max_cells_ = 0;
  }
}

exporter::OrderParaExporter_2::~OrderParaExporter_2() {
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(comm_, &my_rank);
  if (my_rank == 0) {
#endif
    fout_.close();
#ifdef USE_MPI
  }
#endif
}

void exporter::OrderParaExporter_2::coarse_grain(int** n_gl, double** svx_gl, double** svy_gl,
                                                 int& nx, int& ny) const {
  if (nx % 2 != 0 || ny % 2 != 0) {
    std::cout << "Error when coarse grain with nx = " << nx << ", ny = " << ny << std::endl;
    exit(2);
  }
  int nx_new = nx / 2;
  int ny_new = ny / 2;
  int* n_gl_new = new int[nx_new * ny_new]{};
  double* svx_gl_new = new double[nx_new * ny_new]{};
  double* svy_gl_new = new double[nx_new * ny_new]{};
  for (int j = 0; j < ny_new; j++) {
    int j_nx_new = j * nx_new;
    int iy_0_nx = j * 2 * nx;
    int iy_1_nx = iy_0_nx + nx;
    for (int i = 0; i < nx_new; i++) {
      int idx_new = i + j_nx_new;
      int ix_0 = i * 2;
      int ix_1 = ix_0 + 1;
      int idx0 = ix_0 + iy_0_nx;
      int idx1 = ix_1 + iy_0_nx;
      int idx2 = ix_0 + iy_1_nx;
      int idx3 = ix_1 + iy_1_nx;
      n_gl_new[idx_new] = (*n_gl)[idx0] + (*n_gl)[idx1] + (*n_gl)[idx2] + (*n_gl)[idx3];
      svx_gl_new[idx_new] = (*svx_gl)[idx0] + (*svx_gl)[idx1] + (*svx_gl)[idx2] + (*svx_gl)[idx3];
      svy_gl_new[idx_new] = (*svy_gl)[idx0] + (*svy_gl)[idx1] + (*svy_gl)[idx2] + (*svy_gl)[idx3];
    }
  }
  delete[](*n_gl);
  delete[](*svx_gl);
  delete[](*svy_gl);
  *n_gl = n_gl_new;
  *svx_gl = svx_gl_new;
  *svy_gl = svy_gl_new;
  nx = nx_new;
  ny = ny_new;
}

double exporter::OrderParaExporter_2::get_mean_phi(int size, const int* n_gl,
  const double* svx_gl, const double* svy_gl, bool normed) const {
  double phi_sum = 0;
  int count = 0;
  for (int i = 0; i < size; i++) {
    double vx_m = svx_gl[i] / n_gl[i];
    double vy_m = svy_gl[i] / n_gl[i];
    double phi = std::sqrt(vx_m * vx_m + vy_m * vy_m);
    if (normed) {
      phi_sum += phi * n_gl[i];
      count += n_gl[i];
    } else {
      phi_sum += phi;
    }
   
  }
  double phi_mean;
  if (normed) {
    phi_mean = phi_sum / count;
  } else {
    phi_mean = phi_sum / size;
  }
  return phi_mean;
}

void exporter::OrderParaExporter_2::cal_mean_phi(int size, const int* n_gl, const double* svx_gl, const double* svy_gl, double& phi1, double& phi2) {
  phi1 = phi2 = 0;
  int count = 0;
  for (int i = 0; i < size; i++) {
    double vx_m = svx_gl[i] / n_gl[i];
    double vy_m = svy_gl[i] / n_gl[i];
    double phi = std::sqrt(vx_m * vx_m + vy_m * vy_m);
    phi1 += phi * n_gl[i];
    count += n_gl[i];
    phi2 += phi;
  }
  phi1 /= count;
  phi2 /= size;

}

exporter::RhoxExporter::~RhoxExporter() {
  delete[] buf_;
#ifdef USE_MPI
  MPI_File_close(&fh_);
#else
  fout_.close();
#endif
}
