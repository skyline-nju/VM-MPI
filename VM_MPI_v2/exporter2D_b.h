#pragma once
#include "exporter2D.h"
#include "domain2D.h"

class SnapExporter : public BaseExporter {
public:
  explicit SnapExporter(int frame_interval, int n_step);

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

private:
  int count_ = 0;
  char file_prefix_[100];
  MPI_File fh_{};
  int my_rank_;
  int tot_proc_;
};

template<typename TPar>
void SnapExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    char filename[100];
    snprintf(filename, 100, "%s.%04d.bin", file_prefix_, count_);
    count_++;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                  MPI_INFO_NULL, &fh_);
    int my_n = p_arr.size();
    int my_origin;
    int* origin_arr = new int[tot_proc_];
    int* n_arr = new int[tot_proc_];
    MPI_Gather(&my_n, 1, MPI_INT, n_arr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (my_rank_ == 0) {
      origin_arr[0] = 0;
      for (int i = 1; i < tot_proc_; i++) {
        origin_arr[i] = origin_arr[i - 1] + n_arr[i - 1];
      }
    }
    MPI_Scatter(origin_arr, 1, MPI_INT, &my_origin, 1, MPI_INT, 0, MPI_COMM_WORLD);
    delete[] n_arr;
    delete[] origin_arr;

    MPI_Offset offset = my_origin * 3 * sizeof(float);
    float* buf = new float[3 * my_n];
    for (int j = 0; j < my_n; j++) {
      buf[j * 3 + 0] = p_arr[j].pos.x;
      buf[j * 3 + 1] = p_arr[j].pos.y;
      buf[j * 3 + 2] = p_arr[j].theta();
    }
    MPI_File_write_at(fh_, offset, buf, 3 * my_n, MPI_FLOAT, MPI_STATUSES_IGNORE);
    MPI_File_close(&fh_);
    delete[] buf;
  }
}

class RhoxExpoerter : public BaseExporter{
public:
  RhoxExpoerter(int sep, int n_step, const Grid_2& grid, const Vec_2<double>& origin);

  ~RhoxExpoerter() { delete[] buf_; MPI_File_close(&fh_); }

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);
private:
  Vec_2<double> origin_;
  int count_ = 0;
  MPI_File fh_{};
  int my_rank_;
  int tot_proc_;
  MPI_Offset offset_;
  MPI_Offset frame_size_;
  unsigned short* buf_;
  int buf_size_;
};

template<typename TPar>
void RhoxExpoerter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    for (int j = 0; j < buf_size_; j++) {
      buf_[j] = 0;
    }
    const auto end = p_arr.cend();
    for (auto it = p_arr.cbegin(); it != end; ++it) {
      int idx = int((*it).pos.x - origin_.x);
      buf_[idx]++;
    }
    MPI_Offset my_offset = offset_ + frame_size_ * count_;
    MPI_File_write_at(fh_, my_offset, buf_, buf_size_, MPI_UNSIGNED_SHORT, MPI_STATUSES_IGNORE);
    count_++;
  }
}
