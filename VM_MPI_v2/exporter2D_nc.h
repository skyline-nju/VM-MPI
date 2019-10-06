#pragma once

#include "exporter2D.h"
#include "mpi.h"
#include "netcdf.h"
#ifndef _MSC_VER
#include "netcdf_par.h"
#endif

void exporter_nc_ini(int gl_np, double eta0, double eps0, int steps, unsigned long long sd,
                     const Vec_2<double>& gl_l0, const Vec_2<int>& domain_sizes0);

void exporter_nc_finish();

// check whether there is error when outputting netcdf file
void check_err(const int stat, const int line, const char* file);

class NCSnapExporter : public BaseExporter {
public:
  explicit NCSnapExporter(int frame_interval);

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

private:
  char file_prefix_[100];
};

template <typename TPar>
void NCSnapExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    int n = p_arr.size();
    char filename[100];
    snprintf(filename, 100, "%s_%06d.bin", file_prefix_, i_step);
    std::ofstream fout(filename, std::ios::binary);
    double* buf = new double[n * 4];
    for (int i = 0; i < n; i++) {
      buf[i * 4 + 0] = p_arr[i].pos.x;
      buf[i * 4 + 1] = p_arr[i].pos.y;
      buf[i * 4 + 2] = p_arr[i].ori.x;
      buf[i * 4 + 3] = p_arr[i].ori.y;
    }
    fout.write((char*)buf, sizeof(double) * n * 4);
    fout.close();
    delete[] buf;
  }
}

class NCFieldExporter : public BaseExporter {
public:
  explicit NCFieldExporter(int frame_interval,
    const Vec_2<int>& bin_len,
    const Vec_2<int>& domain_rank,
    const Vec_2<int>& gl_cells_size,
    const Vec_2<int>& my_cells_size);

  void set_coarse_grain_box(const Vec_2<int>& gl_cells_size,
    const Vec_2<int>& my_cells_size,
    const Vec_2<int>& domain_rank);

  template <typename TPar, typename T1, typename T2>
  void coarse_grain(const std::vector<TPar>& p_arr, T1* den_fields, T2* vel_fields) const;

  template <typename TPar, typename T1>
  void coarse_grain(const std::vector<TPar>& p_arr, T1* den_fields) const;

  template<typename TPar>
  int dump(int i_step, const std::vector<TPar>& p_arr);

  void open();

private:
  int ncid_;
  int time_id_;
  int spatial_id_;
  int densities_id_;
  int velocities_id_;

  int my_proc_;
  char filename_[100];

  size_t frame_len_;
  size_t spatial_len_ = 2;
  Vec_2<size_t> cg_box_len_;
  size_t gl_field_len_[2];

  size_t den_start_set_[3];
  size_t den_count_set_[3];
  size_t vel_start_set_[4];
  size_t vel_count_set_[4];

  Vec_2<double> origin_;
  size_t time_idx_[1];

  Vec_2<int> n_host_{};

  int n_steps_per_frame_;
};



template <typename TPar, typename T1, typename T2>
void NCFieldExporter::coarse_grain(const std::vector<TPar>& p_arr,
  T1* den_fields, T2* vel_fields) const {
  auto end = p_arr.cend();
  int nx = den_count_set_[2];
  int nx_ny = nx * den_count_set_[1];
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    int ix = int(((*it).pos.x - origin_.x) / cg_box_len_.x);
    int iy = int(((*it).pos.y - origin_.y) / cg_box_len_.y);
    int idx = ix + iy * nx;
    den_fields[idx] += 1;
#ifdef POLAR_ALIGN
    vel_fields[idx] += (*it).ori.x;
    vel_fields[idx + nx_ny] += (*it).ori.y;
#else
    vel_fields[idx] += (*it).ori.x * (*it).ori.x - (*it).ori.y * (*it).ori.y;
    vel_fields[idx + nx_ny] += 2 * (*it).ori.x * (*it).ori.y;
#endif
  }
}

template <typename TPar, typename T1>
void NCFieldExporter::coarse_grain(const std::vector<TPar>& p_arr, T1* den_fields) const {
  auto end = p_arr.cend();
  int nx = den_count_set_[2];
  int nx_ny = nx * den_count_set_[1];
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    int ix = int(((*it).pos.x - origin_.x) / cg_box_len_.x);
    int iy = int(((*it).pos.y - origin_.y) / cg_box_len_.y);
    int idx = ix + iy * nx;
    den_fields[idx] += 1;
  }
}

template <typename TPar>
int NCFieldExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (i_step % n_steps_per_frame_ == 0) {
    open();
    /* inquire varid */
    int stat = nc_inq_varid(ncid_, "time", &time_id_);
    check_err(stat, __LINE__, __FILE__);

    stat = nc_inq_varid(ncid_, "density_field", &densities_id_);
    check_err(stat, __LINE__, __FILE__);

#ifdef OUTPUT_V_FILED
    stat = nc_inq_varid(ncid_, "velocity_field", &velocities_id_);
    check_err(stat, __LINE__, __FILE__);
#endif

#ifndef _MSC_VER
    stat = nc_var_par_access(ncid_, densities_id_, NC_COLLECTIVE);
    check_err(stat, __LINE__, __FILE__);
#ifdef OUTPUT_V_FIELD
    stat = nc_var_par_access(ncid_, velocities_id_, NC_COLLECTIVE);
    check_err(stat, __LINE__, __FILE__);
#endif
    stat = nc_var_par_access(ncid_, time_id_, NC_COLLECTIVE);
    check_err(stat, __LINE__, __FILE__);
#endif
    /* dump variables */
    stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
    check_err(stat, __LINE__, __FILE__);
    time_idx_[0]++;

    const size_t field_size = den_count_set_[1] * den_count_set_[2];
    unsigned short* den_fields = new unsigned short[field_size] {};
#ifdef OUTPUT_V_FILED
    float* vel_fields = new float[2 * field_size]{};
    coarse_grain(p_arr, den_fields, vel_fields);
#else
    coarse_grain(p_arr, den_fields);
#endif

    stat = nc_put_vara(ncid_, densities_id_, den_start_set_, den_count_set_, den_fields);
    check_err(stat, __LINE__, __FILE__);
    den_start_set_[0]++;
    delete[] den_fields;

#ifdef OUTPUT_V_FILED
    stat = nc_put_vara(ncid_, velocities_id_, vel_start_set_, vel_count_set_, vel_fields);
    check_err(stat, __LINE__, __FILE__);
    vel_start_set_[0]++;
    delete[] vel_fields;
#endif
    stat = nc_close(ncid_);
    check_err(stat, __LINE__, __FILE__);
    return 1;
  }
  return 0;
}
