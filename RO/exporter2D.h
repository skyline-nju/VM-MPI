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
#include <map>
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
 * Output the beginning and ending time of the simulation.
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
 *
 * If set use_sub_boxes=1, then calculate local orderparameters over boxes with
 * varies lienar size, otherwise just calcalute globa order parameters.
 */
class OrderParaExporter_2 : public ExporterBase {
public:
  OrderParaExporter_2(const std::string& outfile, int start, int n_step, int sep, 
    const Vec_2<double>& gl_l, MPI_Comm group_comm, int use_sub_boxes=0);

  ~OrderParaExporter_2();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr, int gl_np);

  template <typename TPar>
  void coarse_grain(int** n_gl, double** svx_gl, double** svy_gl, const std::vector<TPar>& p_arr) const;

  void coarse_grain(int** n_gl, double** svx_gl, double** svy_gl, int& nx, int& ny) const;

  double get_mean_phi(int size, const int* n_gl,
    const double* svx_gl, const double* svy_gl, bool normed) const;

  void cal_mean_phi(int size, const int* n_gl, const double* svx_gl, const double* svy_gl,
    double& phi1, double& phi2);

private:
  std::ofstream fout_;
  MPI_Comm comm_;

  std::vector<int> L_arr_;
  bool flag_phi_box_;
  int max_cells_;
};

template<typename TPar>
void OrderParaExporter_2::dump(int i_step, const std::vector<TPar>& p_arr, int gl_np) {
  if (need_export(i_step)) {
    Vec_2<double> gl_vm;
#ifdef USE_MPI
    get_mean_vel(&gl_vm.x, p_arr, gl_np, false, comm_);
    int my_rank, tot_proc;
    MPI_Comm_rank(comm_, &my_rank);
    MPI_Comm_size(comm_, &tot_proc);

    if (flag_phi_box_) {
      int* n_gl = nullptr;
      double* svx_gl = nullptr;
      double* svy_gl = nullptr;
      coarse_grain(&n_gl, &svx_gl, &svy_gl, p_arr);
      int nx = L_arr_.back() / L_arr_[0];
      int ny = nx;

      double phi1, phi2;
      if (my_rank == 0) {
        fout_ << std::fixed << std::setw(10) << std::setprecision(8) << gl_vm.module() << "\t" << atan2(gl_vm.y, gl_vm.x) << "\t";
        cal_mean_phi(nx * ny, n_gl, svx_gl, svy_gl, phi1, phi2);
        fout_ << phi1 << "\t" << phi2;
        while (nx >= 2) {
          coarse_grain(&n_gl, &svx_gl, &svy_gl, nx, ny);
          cal_mean_phi(nx * ny, n_gl, svx_gl, svy_gl, phi1, phi2);
          if (nx == 1) {
            fout_ << "\t" << phi1;
          } else {
            fout_ << "\t" << phi1 << "\t" << phi2;
          }
        }
        fout_ << std::endl;
      }
      delete[] n_gl;
      delete[] svx_gl;
      delete[] svy_gl;
    } else {
      if (my_rank == 0) {
        fout_ << std::fixed << std::setw(16) << std::setprecision(10) << gl_vm.module() << "\t" << atan2(gl_vm.y, gl_vm.x) << std::endl;
      }
    }
#else
    get_mean_vel(&gl_vm.x, p_arr);
    fout_ << gl_vm.module() << "\t" << atan2(gl_vm.y, gl_vm.x) << std::endl;
#endif
  }
}

template<typename TPar>
void OrderParaExporter_2::coarse_grain(int** n_gl, double** svx_gl, double** svy_gl,
  const std::vector<TPar>& p_arr) const {
  int my_rank, tot_proc;
  MPI_Comm_rank(comm_, &my_rank);
  MPI_Comm_size(comm_, &tot_proc);
  //std::cout << "rank = " << my_rank  << "size = " << tot_proc << std::endl;

  int* idx_arr = new int[max_cells_] {};
  int* n_arr = new int[max_cells_] {};
  double* svx_arr = new double[max_cells_] {};
  double* svy_arr = new double[max_cells_] {};
  int* idx_recv{};
  int* n_recv{};
  double* svx_recv{};
  double* svy_recv{};

  double inverse_l = 1. / L_arr_[0];
  int nx = L_arr_.back() / L_arr_[0];

  std::map<int, int> index_map;
  for (const auto& p : p_arr) {
    int idx = int(p.pos.x * inverse_l) + nx * int(p.pos.y * inverse_l);
    auto search = index_map.find(idx);
    int pos;
    if (search != index_map.end()) {
      pos = search->second;
    } else {
      pos = index_map.size();
      index_map.emplace(idx, pos);
    }
    idx_arr[pos] = idx;
    n_arr[pos] += 1;
    svx_arr[pos] += p.ori.x;
    svy_arr[pos] += p.ori.y;
  }
  int size = index_map.size();
  int* size_arr = new int[tot_proc];
  MPI_Gather(&size, 1, MPI_INT, size_arr, 1, MPI_INT, 0, comm_);
  MPI_Bcast(size_arr, tot_proc, MPI_INT, 0, comm_);
  int* displs = new int[tot_proc] {};
  for (int i = 1; i < tot_proc; i++) {
    displs[i] = displs[i - 1] + size_arr[i - 1];
  }
  int gl_size = displs[tot_proc - 1] + size_arr[tot_proc - 1];
  if (my_rank == 0) {
    idx_recv = new int[gl_size];
    n_recv = new int[gl_size];
    svx_recv = new double[gl_size];
    svy_recv = new double[gl_size];
  }

  MPI_Gatherv(idx_arr, size, MPI_INT, idx_recv, size_arr, displs, MPI_INT, 0, comm_);
  MPI_Gatherv(n_arr, size, MPI_INT, n_recv, size_arr, displs, MPI_INT, 0, comm_);
  MPI_Gatherv(svx_arr, size, MPI_DOUBLE, svx_recv, size_arr, displs, MPI_DOUBLE, 0, comm_);
  MPI_Gatherv(svy_arr, size, MPI_DOUBLE, svy_recv, size_arr, displs, MPI_DOUBLE, 0, comm_);
  delete[] size_arr;
  delete[] displs;
  delete[] idx_arr;
  delete[] n_arr;
  delete[] svx_arr;
  delete[] svy_arr;

  if (my_rank == 0) {
    *n_gl = new int[nx * nx]{};
    *svx_gl = new double[nx * nx]{};
    *svy_gl = new double[nx * nx]{};

    for (int i = 0; i < gl_size; i++) {
      int j = idx_recv[i];
      (*n_gl)[j] += n_recv[i];
      (*svx_gl)[j] += svx_recv[i];
      (*svy_gl)[j] += svy_recv[i];
    }
  }
  delete[] idx_recv;
  delete[] n_recv;
  delete[] svx_recv;
  delete[] svy_recv;
}

/**
 * @brief Output snapshot as binary format.
 *
 * For each frame, the information of particles is saved as 3 * N float numbers.
 * 3 float number (x, y, theta) per particle.
 */
class SnapExporter : public ExporterBase {
public:
  explicit SnapExporter(const std::string outfile, int start, int n_step, int sep, int t_beg, MPI_Comm group_comm)
    : ExporterBase(start, n_step, sep), file_prefix_(outfile), t_beg_(t_beg), comm_(group_comm) {}

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

private:
  int count_ = 0;
  std::string file_prefix_;
  int t_beg_;
#ifdef USE_MPI
  MPI_File fh_{};
  MPI_Comm comm_;
#else
  std::ofstream fout_;
#endif
};

template<typename TPar>
void SnapExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    char filename[100];
    snprintf(filename, 100, "%s_%08d.bin", file_prefix_.c_str(), t_beg_ + i_step);
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

/**
 * @brief Exporter for coarse-grained density and momentum field in 2D.
 *
 * Output coarse-grained density and momentum field over boxes with linear
 * size "box_size" every "sep" time steps.
 * 
 * Caution: the number of grids in both x and y directions must be multiple
 * of linear size of boxes for coarse graining.
 */
class FeildExporter : public ExporterBase {
public:
  template <typename TDomain>
  FeildExporter(const std::string& outfile, int start, int n_step, int sep,
    const Grid_2& grid, const TDomain& dm, int box_size = 4);

  template <typename TPar, typename T>
  void coarse_grain(const std::vector<TPar>& p_arr, T* rho, T* vx, T* vy);

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

  void write_data(const float* rho, const float* vx, const float* vy);

  ~FeildExporter();

protected:
  Vec_2<double> origin_;
  Vec_2<int> gl_n_;
  Vec_2<int> n_;
  Vec_2<int> o_;
  Vec_2<double> inv_lc_;

  MPI_File fh_{};
  MPI_Offset frame_size_;
  int n_grids_;
  MPI_Offset* offset_;
  int idx_frame_;
};

template <typename TDomain>
FeildExporter::FeildExporter(const std::string& outfile, int start, int n_step, int sep,
                             const Grid_2& grid, const TDomain& dm, int box_size)
  : ExporterBase(start, n_step, sep), origin_(dm.origin()) {
  gl_n_.x = grid.gl_n().x / box_size;
  gl_n_.y = grid.gl_n().y / box_size;
  n_.x = grid.n().x / box_size;
  n_.y = grid.n().y / box_size;
  o_.x = grid.origin().x / box_size;
  o_.y = grid.origin().y / box_size;
  inv_lc_.x = 1. / box_size;
  inv_lc_.y = 1. / box_size;
  frame_size_ = gl_n_.x * gl_n_.y * sizeof(float) * 3;
  n_grids_ = n_.x * n_.y;

  offset_ = new MPI_Offset[n_.y * 3]{};

  for (int i_field = 0; i_field < 3; i_field++) {
    MPI_Offset start = gl_n_.x * gl_n_.y * sizeof(float) * i_field;
    for (int row = 0; row < n_.y; row++) {
      int idx = row + i_field * n_.y;
      offset_[idx] = start + (o_.x + (o_.y + row) * gl_n_.x) * sizeof(float);
    }
  }

  MPI_File_open(dm.comm(), outfile.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh_);
  MPI_File_set_size(fh_, 0);
  idx_frame_ = 0;
}

template<typename TPar, typename T>
void FeildExporter::coarse_grain(const std::vector<TPar>& p_arr, T* rho, T* vx, T* vy) {
  const auto end = p_arr.end();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    const TPar& p = *it;
    int idx = int((p.pos.x - origin_.x) * inv_lc_.x) + int((p.pos.y - origin_.y) * inv_lc_.y) * n_.x;
    rho[idx] += 1.;
    vx[idx] += p.ori.x;
    vy[idx] += p.ori.y;
  }
}

template<typename TPar>
void FeildExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    float* rho = new float[n_grids_]();
    float* vx = new float[n_grids_]();
    float* vy = new float[n_grids_]();
    coarse_grain(p_arr, rho, vx, vy);

    write_data(rho, vx, vy);
    
    delete[] rho;
    delete[] vx;
    delete[] vy;
  }
}

/**
 * @brief Output time-averaged fields of density and momentum.
 * 
 * Accumulate density and momentum fields over "sep" time steps, then cal their
 * time aveage.
 * Caution: the number of grids in both x and y directions must be multiple of
 * linear size of boxes for coarse graining.
 */
class TimeAveFeildExporter : public FeildExporter {
public:
  template <typename TDomain>
  TimeAveFeildExporter(const std::string& outfile, int start, int n_step, int sep,
    const Grid_2& grid, const TDomain& dm, int box_size = 4);

  ~TimeAveFeildExporter();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);
  
private:
  int t_win_;
  double* sum_n_;
  double* sum_vx_;
  double* sum_vy_;
};

template<typename TDomain>
TimeAveFeildExporter::TimeAveFeildExporter(const std::string& outfile, int start, int n_step, int sep,
                                           const Grid_2& grid, const TDomain& dm, int box_size)
  : FeildExporter(outfile, start, n_step, sep, grid, dm, box_size), t_win_(sep){
  sum_n_ = new double[n_grids_]();
  sum_vx_ = new double[n_grids_]();
  sum_vy_ = new double[n_grids_]();
}

template<typename TPar>
void TimeAveFeildExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (i_step > start_) {
    coarse_grain(p_arr, sum_n_, sum_vx_, sum_vy_);
    if ((i_step - start_) % t_win_ == 0 && i_step <= start_ + n_step_) {
      float* rho = new float[n_grids_];
      float* vx = new float[n_grids_];
      float* vy = new float[n_grids_];
      for (int i = 0; i < n_grids_; i++) {
        rho[i] = sum_n_[i] / t_win_;
        vx[i] = sum_vx_[i] / t_win_;
        vy[i] = sum_vy_[i] / t_win_;
      }
      write_data(rho, vx, vy);
      delete[] rho;
      delete[] vx;
      delete[] vy;
      for (int i = 0; i < n_grids_; i++) {
        sum_n_[i] = sum_vx_[i] = sum_vy_[i] = 0.;
      }
    }
  }
}

}
