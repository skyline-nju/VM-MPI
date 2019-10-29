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

exporter::LogExporter::LogExporter(const std::string& outfile, int start, int n_step, int sep, int np)
  : ExporterBase(start, n_step, sep), n_par_(np) {
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
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
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
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
    fout << "speed=" << std::scientific << n_step_ * double(n_par_) / elapsed_seconds.count()
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
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
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
}

exporter::OrderParaExporter_2::OrderParaExporter_2(const std::string& outfile, int start, int n_step, int sep)
  : ExporterBase(start, n_step, sep) {
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0) {
#endif
    fout_.open(outfile);
#ifdef USE_MPI
  }
#endif

}

exporter::OrderParaExporter_2::~OrderParaExporter_2() {
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0) {
#endif
    fout_.close();
#ifdef USE_MPI
  }
#endif
}

exporter::RhoxExporter::~RhoxExporter() {
  delete[] buf_;
#ifdef USE_MPI
  MPI_File_close(&fh_);
#else
  fout_.close();
#endif
}
