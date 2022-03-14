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
#include "node.h"

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
  //ExporterBase() : n_step_(0) {}

  ExporterBase(int start, int n_step, int sep, MPI_Comm group_comm);

  void set_lin_frame(int start, int n_step, int sep);

  bool need_export(const int i_step);

protected:
  int n_step_;    // total steps to run
  int start_ = 0; // The first step
  int my_rank_=0;
  int tot_proc_=1;
#ifdef USE_MPI
  MPI_Comm comm_;
#endif
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
  LogExporter(const std::string& outfile,
              int start, int n_step, int sep,
              int np, 
              MPI_Comm group_comm);
  ~LogExporter();

  void record(int i_step);

  std::ofstream fout;
private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  int n_par_;
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
  OrderParaExporter_2(const std::string& outfile, 
                      int start, int n_step, int sep, 
                      const Vec_2<double>& gl_l,
                      MPI_Comm group_comm);

  ~OrderParaExporter_2();

  void dump(int i_step, const std::vector<BiNode<Bird_2>>& p_arr, int gl_np);

private:
  std::ofstream fout_;
};


/**
 * @brief Output snapshot as binary format.
 *
 * For each frame, the information of particles is saved as 3 * N float numbers.
 * 3 float number (x, y, theta) per particle.
 */
class SnapExporter : public ExporterBase {
public:
  explicit SnapExporter(const std::string outfile,
                        int start, int n_step, int sep,
                        int t_beg, 
                        MPI_Comm group_comm)
    : ExporterBase(start, n_step, sep, group_comm), 
      file_prefix_(outfile), t_beg_(t_beg) {}

  void dump(int i_step, const std::vector<BiNode<Bird_2>>& p_arr);

private:
  int count_ = 0;
  std::string file_prefix_;
  int t_beg_;
#ifdef USE_MPI
  MPI_File fh_{};
#else
  std::ofstream fout_;
#endif
};


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
  int n_grids_;

#ifdef USE_MPI
  MPI_File fh_{};
  MPI_Offset frame_size_;
  MPI_Offset* offset_;
#else
  std::ofstream fout_;
#endif
  int idx_frame_;
};

template <typename TDomain>
FeildExporter::FeildExporter(const std::string& outfile, int start, int n_step, int sep,
                             const Grid_2& grid, const TDomain& dm, int box_size)
  : ExporterBase(start, n_step, sep, dm.comm()), origin_(dm.origin()) {
  gl_n_.x = grid.gl_n().x / box_size;
  gl_n_.y = grid.gl_n().y / box_size;
  n_.x = grid.n().x / box_size;
  n_.y = grid.n().y / box_size;
  o_.x = grid.origin().x / box_size;
  o_.y = grid.origin().y / box_size;
  inv_lc_.x = 1. / box_size;
  inv_lc_.y = 1. / box_size;
  n_grids_ = n_.x * n_.y;

#ifdef USE_MPI
  frame_size_ = gl_n_.x * gl_n_.y * sizeof(float) * 3;
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
#else
  fout_.open(outfile.c_str(), std::ios::binary);

#endif
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

void create_folders(int my_rank, MPI_Comm group_comm);

}
