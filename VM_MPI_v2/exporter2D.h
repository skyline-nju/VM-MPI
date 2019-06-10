#pragma once

#include "config.h"
#include "comn.h"
#include "vect.h"
#include "netcdf.h"
#ifndef _MSC_VER
#include "netcdf_par.h"
#endif
#include "mpi.h"
#include "cellList2D.h"

void ini_output(int gl_np, double eta0, double eps0, int steps, unsigned long long sd,
                const Vec_2<double> &gl_l0, const Vec_2<int> &domain_sizes0);

void output_finalize();

// check whether there is error when outputting netcdf file
void check_err(const int stat, const int line, const char * file);

int get_start_particle_num(int particle_num);

template <typename TPar>
void get_vel_sum(double *vel_sum, const std::vector<TPar> &p_arr) {
  vel_sum[0] = vel_sum[1] = 0.;
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it!=end; ++it) {
#ifdef POLAR_ALIGN
    vel_sum[0] += (*it).ori.x;
    vel_sum[1] += (*it).ori.y;
#else
    vel_sum[0] += (*it).ori.x * (*it).ori.x - (*it).ori.y * (*it).ori.y;
    vel_sum[1] += 2 * (*it).ori.x * (*it).ori.y;
#endif
  }
}

template <typename TPar>
void get_mean_vel(double *vel_mean, const std::vector<TPar> &p_arr,
                  int gl_np, bool flag_broadcast) {
  get_vel_sum(vel_mean, p_arr);
  double gl_vel_sum[2];
  MPI_Reduce(vel_mean, gl_vel_sum, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0) {
    vel_mean[0] = gl_vel_sum[0] / gl_np;
    vel_mean[1] = gl_vel_sum[1] / gl_np;
  }
  if (flag_broadcast) {
    MPI_Bcast(vel_mean, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
}

class LogExporter : public BaseLogExporter {
public:
  explicit LogExporter(int interval);
};

class OrderParaExporter : public BaseExporter {
public:
  explicit OrderParaExporter(int interval);

  ~OrderParaExporter();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar> &p_arr);

  template <typename TPar, typename T>
  void dump(int i_step, const std::vector<TPar> &p_arr, CellListNode_2<TPar> &cl,
            const std::vector<T> &n_arr);
private:
  std::ofstream fout_;
  int gl_np_;
  int my_proc_;
};

template <typename TPar>
void OrderParaExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    Vec_2<double> gl_vm;
    get_mean_vel(&gl_vm.x, p_arr, gl_np_, false);
    if (my_proc_ == 0) {
      fout_ << gl_vm.module() << "\t" << gl_vm.x << "\t" << gl_vm.y << std::endl;
    }
  }
}

template <typename TPar, typename T>
void OrderParaExporter::dump(int i_step, const std::vector<TPar>& p_arr,
                             CellListNode_2<TPar>& cl, const std::vector<T>& n_arr) {
  if (need_export(i_step)) {
    Vec_2<double> gl_vm{};
    int np = cl.get_par_num(n_arr);
    int gl_np;
    MPI_Reduce(&np, &gl_np, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    get_mean_vel(&gl_vm.x, p_arr, gl_np, false);
    if (my_proc_ == 0) {
      fout_ << gl_vm.module() << "\t" << gl_vm.x << "\t" << gl_vm.y << "\t" << gl_np;
      if (i_step % 5000 == 0) {
        fout_ << std::endl;
      } else {
        fout_ << "\n";
      }
    }
  }
}

class FieldExporter : public BaseExporter {
public:
  explicit FieldExporter(int frame_interval,
                         int bin_len,
                         const Vec_2<int> &domain_rank,
                         const Vec_2<int> &gl_cells_size,
                         const Vec_2<int> &my_cells_size);

  void set_coarse_grain_box(const Vec_2<int> &gl_cells_size,
                            const Vec_2<int> &my_cells_size,
                            const Vec_2<int> &domain_rank);

  template <typename TPar, typename T1, typename T2>
  void coarse_grain(const std::vector<TPar> &p_arr, T1 *den_fields, T2 *vel_fields) const;

  template<typename TPar>
  int dump(int i_step, const std::vector<TPar> &p_arr);

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
  size_t cg_box_len_;
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
void FieldExporter::coarse_grain(const std::vector<TPar>& p_arr,
                                 T1* den_fields, T2* vel_fields) const {
  auto end = p_arr.cend();
  int nx = den_count_set_[2];
  int nx_ny = nx * den_count_set_[1];
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    int ix = int(((*it).pos.x - origin_.x) / cg_box_len_);
    int iy = int(((*it).pos.y - origin_.y) / cg_box_len_);
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

template <typename TPar>
int FieldExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (i_step % n_steps_per_frame_ == 0) {
    open();
    /* inquire varid */
    int stat = nc_inq_varid(ncid_, "time", &time_id_);
    check_err(stat, __LINE__, __FILE__);

    stat = nc_inq_varid(ncid_, "density_field", &densities_id_);
    check_err(stat, __LINE__, __FILE__);

    stat = nc_inq_varid(ncid_, "velocity_field", &velocities_id_);
    check_err(stat, __LINE__, __FILE__);

#ifndef _MSC_VER
    stat = nc_var_par_access(ncid_, densities_id_, NC_COLLECTIVE);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_var_par_access(ncid_, velocities_id_, NC_COLLECTIVE);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_var_par_access(ncid_, time_id_, NC_COLLECTIVE);
    check_err(stat, __LINE__, __FILE__);
#endif
    /* dump variables */
    stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
    check_err(stat, __LINE__, __FILE__);

    const size_t field_size = den_count_set_[1] * den_count_set_[2];
    unsigned short* den_fields = new unsigned short[field_size] {};
    float* vel_fields = new float[2 * field_size]{};

    coarse_grain(p_arr, den_fields, vel_fields);

    stat = nc_put_vara(ncid_, densities_id_, den_start_set_, den_count_set_, den_fields);
    check_err(stat, __LINE__, __FILE__);

    stat = nc_put_vara(ncid_, velocities_id_, vel_start_set_, vel_count_set_, vel_fields);
    check_err(stat, __LINE__, __FILE__);
    time_idx_[0]++;
    den_start_set_[0]++;
    vel_start_set_[0]++;

    delete[] den_fields;
    delete[] vel_fields;
    stat = nc_close(ncid_);
    check_err(stat, __LINE__, __FILE__);
    return 1;
  }
  return 0;
}