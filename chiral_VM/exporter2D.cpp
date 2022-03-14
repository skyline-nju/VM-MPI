#include "exporter2D.h"

using namespace exporter;
ExporterBase::ExporterBase(int start, int n_step, int sep, 
                           MPI_Comm group_comm) {
#ifdef USE_MPI
  comm_ = group_comm;
  MPI_Comm_rank(comm_, &my_rank_);
  MPI_Comm_size(comm_, &tot_proc_);
#endif
  set_lin_frame(start, n_step, sep);
}


void ExporterBase::set_lin_frame(int start, int n_step, int sep) {
  n_step_ = n_step;
  start_ = start;
  for (auto i = start + sep; i <= n_step_; i += sep) {
    frames_arr_.push_back(i);
  }
  frame_iter_ = frames_arr_.begin();
}

bool ExporterBase::need_export(int i_step) {
  bool flag = false;
  if (!frames_arr_.empty() && i_step == (*frame_iter_)) {
    frame_iter_++;
    flag = true;
  }
  return flag;
}

LogExporter::LogExporter(const std::string& outfile, 
                         int start, int n_step, int sep,
                         int np, 
                         MPI_Comm group_comm)
  : ExporterBase(start, n_step, sep, group_comm), n_par_(np) {
  if (my_rank_ == 0) {
    fout.open(outfile);
    t_start_ = std::chrono::system_clock::now();
    auto start_time = std::chrono::system_clock::to_time_t(t_start_);
    char str[100];
    std::strftime(str, 100, "%c", std::localtime(&start_time));
    fout << "Started simulation at " << str << "\n";
  }
}

LogExporter::~LogExporter() {
  if (my_rank_ == 0) {
    const auto t_now = std::chrono::system_clock::now();
    auto end_time = std::chrono::system_clock::to_time_t(t_now);
    char str[100];
    // ReSharper disable CppDeprecatedEntity
    std::strftime(str, 100, "%c", std::localtime(&end_time));
    // ReSharper restore CppDeprecatedEntity
    fout << "Finished simulation at " << str << "\n";
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    fout << "speed=" << std::scientific << step_count_ * double(n_par_) / elapsed_seconds.count() /tot_proc_
      << " particle time step per second per core\n";
    fout.close();
  }
}

void LogExporter::record(int i_step) {
  if (need_export(i_step) && my_rank_ == 0) {
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


OrderParaExporter_2::OrderParaExporter_2(const std::string& outfile,
                                         int start, int n_step, int sep,
                                         const Vec_2<double>& gl_l,
                                         MPI_Comm group_comm)
  : ExporterBase(start, n_step, sep, group_comm){
  if (my_rank_ == 0) {
    fout_.open(outfile);
  }
}

OrderParaExporter_2::~OrderParaExporter_2() {
  if (my_rank_ == 0) {
    fout_.close();
  }
}

void OrderParaExporter_2::dump(int i_step,
                               const std::vector<BiNode<Bird_2>>& p_arr,
                               int gl_np) {
  if (need_export(i_step)) {
    double my_v[2]{};
    for (const auto& p : p_arr) {
      my_v[0] += p.ori.x;
      my_v[1] += p.ori.y;
    }
    double gl_v[2]{};
#ifdef USE_MPI
    MPI_Reduce(my_v, gl_v, 2, MPI_DOUBLE, MPI_SUM, 0, comm_);
#else
    gl_v[0] = my_v[0];
    gl_v[1] = my_v[1];
#endif
    if (my_rank_ == 0) {
      gl_v[0] /= gl_np;
      gl_v[1] /= gl_np;
    }
    double m = std::sqrt(gl_v[0] * gl_v[0] + gl_v[1] * gl_v[1]);
    double theta = std::atan2(gl_v[1], gl_v[0]);

    if (my_rank_ == 0) {
      fout_ << std::fixed << std::setw(16) << std::setprecision(10)
        << m << "\t" << theta << "\t" << std::endl;
    }
  }
}

void SnapExporter::dump(int i_step, const std::vector<BiNode<Bird_2>>& p_arr) {
  if (need_export(i_step)) {
    char filename[200];
    snprintf(filename, 200, "%s_%08d.bin", file_prefix_.c_str(), t_beg_ + i_step);
    count_++;
    int my_n = p_arr.size();
#ifdef USE_MPI
    MPI_File_open(comm_, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                  MPI_INFO_NULL, &fh_);
    int my_rank;
    MPI_Comm_rank(comm_, &my_rank);
    int tot_proc;
    MPI_Comm_size(comm_, &tot_proc);
    int my_origin;
    int* origin_arr = new int[tot_proc];
    int* n_arr = new int[tot_proc];
    MPI_Gather(&my_n, 1, MPI_INT, n_arr, 1, MPI_INT, 0, comm_);
    if (my_rank == 0) {
      origin_arr[0] = 0;
      for (int i = 1; i < tot_proc; i++) {
        origin_arr[i] = origin_arr[i - 1] + n_arr[i - 1];
      }
    }
    MPI_Scatter(origin_arr, 1, MPI_INT, &my_origin, 1, MPI_INT, 0, comm_);
    delete[] n_arr;
    delete[] origin_arr;

    MPI_Offset offset = my_origin * 3 * sizeof(float);
#else
    fout_.open(filename, std::ios::binary);
#endif
    float* buf = new float[3 * my_n];
    for (int j = 0; j < my_n; j++) {
      buf[j * 3 + 0] = p_arr[j].pos.x;
      buf[j * 3 + 1] = p_arr[j].pos.y;
      buf[j * 3 + 2] = p_arr[j].get_theta();
    }

#ifdef USE_MPI
    MPI_File_write_at(fh_, offset, buf, 3 * my_n, MPI_FLOAT, MPI_STATUSES_IGNORE);
    MPI_File_close(&fh_);
#else
    fout_.write((char*)buf, sizeof(float) * my_n * 3);
    fout_.close();
#endif
    delete[] buf;
  }
}

FeildExporter::~FeildExporter() {
#ifdef USE_MPI
  delete[] offset_;
  MPI_File_close(&fh_);
#else
  fout_.close();
#endif
}

void FeildExporter::write_data(const float* rho, const float* vx, const float* vy) {
#ifdef USE_MPI
  const MPI_Offset frame_start = idx_frame_ * frame_size_;
  MPI_Offset offset1, offset2, offset3;
  int pos = 0;
  for (int row = 0; row < n_.y; row++) {
    offset1 = frame_start + offset_[row];
    offset2 = frame_start + offset_[row + n_.y];
    offset3 = frame_start + offset_[row + n_.y + n_.y];
    MPI_File_write_at(fh_, offset1, &rho[pos], n_.x, MPI_FLOAT, MPI_STATUSES_IGNORE);
    MPI_File_write_at(fh_, offset2, &vx[pos],  n_.x, MPI_FLOAT, MPI_STATUSES_IGNORE);
    MPI_File_write_at(fh_, offset3, &vy[pos],  n_.x, MPI_FLOAT, MPI_STATUSES_IGNORE);
    pos += n_.x;
  }
#else
  size_t buf_size = sizeof(float) * n_grids_;
  fout_.write((char*)rho, buf_size);
  fout_.write((char*)vx, buf_size);
  fout_.write((char*)vy, buf_size);
#endif
  idx_frame_++;
}

void exporter::create_folders(int my_rank, MPI_Comm group_comm) {
  // output setting
  if (my_rank == 0) {
    mkdir("data");
#ifdef _MSC_VER
    mkdir("data\\snap");
#else
    mkdir("data/snap");
#endif
  }

#ifdef USE_MPI
  MPI_Barrier(group_comm);
#endif
}
