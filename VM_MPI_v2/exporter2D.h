/**
 * @file exporter2D.h
 * @author Yu Duan (duanyu100@yeah.net)
 * @brief output logfile, order parameters, snapshot in binary format and so on.
 * @version 0.1
 * @date 2019-10-29
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once
#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "config.h"
#include "particle2D.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace exporter {

#ifdef _MSC_VER
const std::string delimiter("\\");
#else
const std::string delimiter("/");
#endif

/**
 * @brief Basic class for exporting data.
 * 
 * Define the timming to dump data. 
 */
class ExporterBase {
public:
  ExporterBase() : n_step_(0) {}

  ExporterBase(int start, int n_step, int sep) : start_(start), n_step_(n_step) {
    set_lin_frame(start, n_step, sep);
  }

  void set_lin_frame(int start, int n_step, int sep);

  bool need_export(const int i_step);

protected:
  int n_step_;    // total steps to run
  int start_ = 0; // The first step 
private:
  std::vector<int> frames_arr_; // frames that need to export
  std::vector<int>::iterator frame_iter_;
};

/**
 * @brief Exporter to output log
 * 
 * Output the parameters after the initialization.
 * Output the beginning and endding time of the simulation.
 * Record time every certain time steps.
 */
class LogExporter : public ExporterBase {
public:
  LogExporter(const std::string& outfile, int start, int n_step, int sep, int np, MPI_Comm group_comm);

  ~LogExporter();

  void record(int i_step);

  std::ofstream fout;
private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  int n_par_;
  MPI_Comm comm_;
  int step_count_ = 0;
};

/**
 * @brief Output order parameters.
 */
class OrderParaExporter_2 : public ExporterBase {
public:
  OrderParaExporter_2(const std::string& outfile, int start, int n_step, int sep, MPI_Comm group_comm);

  ~OrderParaExporter_2();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr, int gl_np);
private:
  std::ofstream fout_;
  MPI_Comm comm_;
};

template<typename TPar>
void OrderParaExporter_2::dump(int i_step, const std::vector<TPar>& p_arr, int gl_np) {
  if (need_export(i_step)) {
    Vec_2<double> gl_vm;
#ifdef USE_MPI
    get_mean_vel(&gl_vm.x, p_arr, gl_np, false, comm_);
    int my_rank;
    MPI_Comm_rank(comm_, &my_rank);
    if (my_rank == 0) {
      fout_ << std::fixed << std::setw(16) << std::setprecision(10) << gl_vm.module() << "\t" << atan2(gl_vm.y, gl_vm.x) << std::endl;
    }
#else
    get_mean_vel(&gl_vm.x, p_arr);
    fout_ << gl_vm.module() << "\t" << atan2(gl_vm.y, gl_vm.x) << std::endl;
#endif
  }
}

/**
 * @brief Output snapshot as binary format.
 * 
 * For each frame, the information of particles is saved as 3 * N float numbers.
 * 3 float number (x, y, theta) per particle.
 */
class SnapExporter : public ExporterBase {
public:
  explicit SnapExporter(const std::string outfile, int start, int n_step, int sep, MPI_Comm group_comm)
    : ExporterBase(start, n_step, sep), file_prefix_(outfile), comm_(group_comm) {}

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

private:
  int count_ = 0;
  std::string file_prefix_;
#ifdef USE_MPI
  MPI_File fh_{};
  MPI_Comm comm_;
#else
  ofsteam fout_;
#endif
};

template<typename TPar>
void SnapExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    char filename[100];
    snprintf(filename, 100, "%s.%04d.bin", file_prefix_.c_str(), count_);
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
      buf[j * 3 + 2] = p_arr[j].theta();
    }

#ifdef USE_MPI
    MPI_File_write_at(fh_, offset, buf, 3 * my_n, MPI_FLOAT, MPI_STATUSES_IGNORE);
    MPI_File_close(&fh_);
#else
    fout_.write(buf, sizeof(float) * my_n * 3);
    fout_.close()
#endif
    delete[] buf;
  }
}

/**
 * @brief Output density profile averaged along y direction, rho_y(x)
 */
class RhoxExporter : public ExporterBase {
public:
  template <typename TDomain>
  RhoxExporter(const std::string& outfile, int start, int n_step, int sep,
                const Grid_2& grid, const TDomain& dm);

  ~RhoxExporter();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);
private:
  Vec_2<double> origin_;
  int count_ = 0;
#ifdef USE_MPI
  MPI_File fh_{};
  int my_rank_;
  int tot_proc_;
  MPI_Offset offset_;
  MPI_Offset frame_size_;
#else
  std::fstream fout_;
#endif
  float* buf_;
  int buf_size_;
  double bin_area_;
};

template <typename TDomain>
RhoxExporter::RhoxExporter(const std::string& outfile, int start, int n_step, int sep,
                           const Grid_2& grid, const TDomain& dm)
  : ExporterBase(start, n_step, sep), origin_(dm.origin()) {
  frame_size_ = grid.gl_n().x * sizeof(float);
  buf_size_ = grid.n().x;
  buf_ = new float[buf_size_];
  bin_area_ = dm.gl_l().y;
#ifdef USE_MPI
  MPI_Comm_size(dm.comm(), &tot_proc_);
  MPI_Comm_rank(dm.comm(), &my_rank_);
  offset_ = grid.origin().x * sizeof(float);
  MPI_File_open(dm.comm(), outfile.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE,
    MPI_INFO_NULL, &fh_);
#else
    fout_.open(outfile.c_str(), std::ios::binary);
#endif
}

template<typename TPar>
void RhoxExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    for (int j = 0; j < buf_size_; j++) {
      buf_[j] = 0;
    }
    const auto end = p_arr.cend();
    for (auto it = p_arr.cbegin(); it != end; ++it) {
      int idx = int((*it).pos.x - origin_.x);
      buf_[idx] += 1;
    }

    for (int j = 0; j < buf_size_; j++) {
      buf_[j] /= bin_area_;
    }
#ifdef USE_MPI
    MPI_Offset my_offset = offset_ + frame_size_ * count_;
    MPI_File_write_at(fh_, my_offset, buf_, buf_size_, MPI_FLOAT, MPI_STATUSES_IGNORE);
#else
    fout_.write(buf, sizeof(float) * buf_size_);
#endif
    count_++;
  }
}

}