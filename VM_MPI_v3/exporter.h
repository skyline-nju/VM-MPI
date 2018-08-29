#pragma once

#include "comn.h"
#include "vect.h"
#include "netcdf.h"
#ifndef _MSC_VER
#include "netcdf_par.h"
#endif
#include "mpi.h"

void ini_output(int gl_np, double eta0, double eps0, int steps, unsigned long long sd,
                const Vec_3<double> &gl_l0, const Vec_3<int> &domain_sizes0);

// check whether there is error when outputting netcdf file
void check_err(const int stat, const int line, const char * file);

int get_start_particle_num(int particle_num);

template <typename T>
void get_new_axis(const Vec_3<T> &new_x_axis, Vec_3<T> &new_y_axis, Vec_3<T> &new_z_axis) {
  if (new_x_axis.x == 0. && new_x_axis.y == 0.) {
    new_z_axis = Vec_3<double>(-1, 0, 0);
    new_y_axis = Vec_3<double>(0, 1, 0);
  } else {
    new_y_axis = Vec_3<double>(-new_x_axis.y, new_x_axis.x, 0);
    new_y_axis.normalize();
    new_z_axis = new_x_axis.cross(new_y_axis);
    //double length = new_x_axis.module();
    //double length_xy = std::sqrt(new_x_axis.x * new_x_axis.x + new_x_axis.y * new_x_axis.y);
    //double c = length_xy / length;
    //double s = new_x_axis.z / length;
    //new_z_axis = Vec_3<double>(0, 0, 1);
    //new_z_axis.rotate(c, -s, new_y_axis);
  }
}

template <typename TPar>
void get_mean_vel(double *vel_mean, const std::vector<TPar> p_arr,
                  int gl_np, bool flag_broadcast) {
  vel_mean[0] = vel_mean[1] = vel_mean[2] = 0;
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    vel_mean[0] += (*it).ori.x;
    vel_mean[1] += (*it).ori.y;
    vel_mean[2] += (*it).ori.z;
  }
  double gl_vel_sum[3];
  MPI_Reduce(vel_mean, gl_vel_sum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0) {
    vel_mean[0] = gl_vel_sum[0] / gl_np;
    vel_mean[1] = gl_vel_sum[1] / gl_np;
    vel_mean[2] = gl_vel_sum[2] / gl_np;
  }
  if (flag_broadcast) {
    MPI_Bcast(vel_mean, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
}

class LogExporter: public BaseLogExporter {
public:
  explicit LogExporter(int interval);
};

class OrderParaExporter: public BaseExporter {
public:
  explicit OrderParaExporter(int interval);

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar> &p_arr);
private:
  std::ofstream fout_;
  int gl_np_;
  int my_proc_;
};


template <typename TPar>
void OrderParaExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if(need_export(i_step)) {
    Vec_3<double> gl_vm;
    get_mean_vel(&gl_vm.x, p_arr, gl_np_, false);
    if (my_proc_ == 0) {
      fout_ << gl_vm.module() << "\t" << gl_vm.x << "\t" << gl_vm.y << "\t" << gl_vm.z << std::endl;
    }
  }
}

class FieldExporter: public BaseExporter {
public:
  explicit FieldExporter(int frame_interval, int first_frame, int bin_len,
                        const Vec_3<int> &domain_rank,
                        const Vec_3<int> &gl_cells_size,
                        const Vec_3<int> &my_cells_size);

  void set_coarse_grain_box(const Vec_3<int> &gl_cells_size,
                            const Vec_3<int> &my_cells_size,
                            const Vec_3<int> &domain_rank);

  template <typename TPar, typename T1, typename T2>
  void coarse_grain(const std::vector<TPar> &p_arr, T1 *den_fields, T2 *vel_fields) const;

  template<typename TPar>
  void dump(int i_step, const std::vector<TPar> &p_arr);

private:
  int ncid_;
  int time_id_;
  int spatial_id_;
  int densities_id_;
  int velocities_id_;

  int my_proc_;

  size_t frame_len_;
  size_t spatial_len_ = 3;
  size_t cg_box_len_;
  size_t gl_field_len_[3];

  size_t den_start_set_[4];
  size_t den_count_set_[4];
  size_t vel_start_set_[5];
  size_t vel_count_set_[5];

  Vec_3<double> origin_;
  size_t time_idx_[1];
};


template <typename TPar, typename T1, typename T2>
void FieldExporter::coarse_grain(const std::vector<TPar>& p_arr,
                                T1* den_fields, T2* vel_fields) const {
  auto end = p_arr.cend();
  int nx = den_count_set_[3];
  int nx_ny = nx * den_count_set_[2];
  int nx_ny_nz = nx_ny * den_count_set_[1];
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    int ix = int(((*it).pos.x - origin_.x) / cg_box_len_);
    int iy = int(((*it).pos.y - origin_.y) / cg_box_len_);
    int iz = int(((*it).pos.z - origin_.z) / cg_box_len_);
    int idx = ix + iy * nx + iz * nx_ny;
    den_fields[idx] += 1;
    vel_fields[idx] += (*it).ori.x;
    vel_fields[idx + nx_ny_nz] += (*it).ori.y;
    vel_fields[idx + nx_ny_nz + nx_ny_nz] += (*it).ori.z;
  }
}

template <typename TPar>
void FieldExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    /* time step */
    int stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
    check_err(stat, __LINE__, __FILE__);

    /* fields */
    const size_t field_size = den_count_set_[1] * den_count_set_[2] * den_count_set_[3];
    unsigned short* den_fields = new unsigned short[field_size] {};
    float* vel_fields = new float[3 * field_size]{};

    coarse_grain(p_arr, den_fields, vel_fields);
    
    //std::cout << "density start set\t" << den_start_set_[0] << "\t" << den_start_set_[1] << "\t" << den_start_set_[2] << "\t" << den_start_set_[3] << std::endl;
    //std::cout << "density count set\t" << den_count_set_[0] << "\t" << den_count_set_[1] << "\t" << den_count_set_[2] << "\t" << den_count_set_[3] << std::endl;
    stat = nc_put_vara(ncid_, densities_id_, den_start_set_, den_count_set_, den_fields);
    check_err(stat, __LINE__, __FILE__);

    //std::cout << "velocity start set\t" << vel_start_set_[0] << "\t" << vel_start_set_[1] << "\t" << vel_start_set_[2] << "\t" << vel_start_set_[3] << "\t" << vel_start_set_[4] << std::endl;
    //std::cout << "velocity count set\t" << vel_count_set_[0] << "\t" << vel_count_set_[1] << "\t" << vel_count_set_[2] << "\t" << vel_count_set_[3] << "\t" << vel_count_set_[4] << std::endl;
    stat = nc_put_vara(ncid_, velocities_id_, vel_start_set_, vel_count_set_, vel_fields);
    check_err(stat, __LINE__, __FILE__);
    time_idx_[0]++;
    den_start_set_[0]++;
    vel_start_set_[0]++;

    nc_sync(ncid_);
    delete[] den_fields;
    delete[] vel_fields;
  }
}

class ParticleExporter: public BaseExporter {
public:
  explicit ParticleExporter(int frame_interval, int first_frame,
                            bool flag_vel = false, bool flag_ori=true);

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar> &p_arr);

private:
  /* id for each variables */
  int ncid_;
  int time_id_;
  int spatial_id_;
  int cell_spatial_id_;
  int cell_lengths_id_;
  int cell_origin_id_;
  int coordinates_id_;
  int velocities_id_;

  int vel_mean_id_;
  int vel_parallel_id_;
  int theta_id_;

  /* dimension lengths */
  const size_t frame_len_;
  const size_t spatial_len_ = 3;
  const size_t atom_len_;
  const size_t cell_spatial_len_ = 3;
  float cell_origin_data_[3];
  float cell_lengths_data_[3];
  size_t time_idx_[1];

  bool flag_vel_;
  bool flag_ori_;
  int gl_np_;
};

template<typename TPar>
void ParticleExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    /* time step */
    int stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
    check_err(stat, __LINE__, __FILE__);

    /* cell origin and lengths */
    {
      size_t startset[2] = { time_idx_[0], 0 };
      size_t countset[2] = { 1, cell_spatial_len_ };
      stat = nc_put_vara(ncid_, cell_lengths_id_, startset, countset, cell_lengths_data_);
      check_err(stat, __LINE__, __FILE__);
      stat = nc_put_vara(ncid_, cell_origin_id_, startset, countset, cell_origin_data_);
      check_err(stat, __LINE__, __FILE__);
    }
    /* coordinates */
    {
      int particle_count = p_arr.size();
      int particle_begin = get_start_particle_num(particle_count);

      size_t startset[3] = { time_idx_[0], particle_begin, 0 };
      size_t countset[3] = { 1, particle_count, spatial_len_ };
      
      float *coordinates_buf = new float[spatial_len_ * particle_count];
      float *velocities_buf = nullptr;
      Vec_3<double> vel_mean{};
      Vec_3<double> z_axis_new(0, 0, 1);
      Vec_3<double> y_axis_new(0, 1, 0);
      Vec_3<double> x_axis_new;
      float * vel_para_module = nullptr;
      float * theta = nullptr;
      if (flag_vel_) {
        velocities_buf = new float[spatial_len_ * particle_count];
      }
      if (flag_ori_) {
        vel_para_module = new float[particle_count];
        theta = new float[particle_count];
        get_mean_vel(&vel_mean.x, p_arr, gl_np_, true);
        x_axis_new = vel_mean / vel_mean.module();
        get_new_axis(x_axis_new, y_axis_new, z_axis_new);
      }
      for (int i = 0; i < particle_count; i++) {
        coordinates_buf[i * 3    ] = p_arr[i].pos.x;
        coordinates_buf[i * 3 + 1] = p_arr[i].pos.y;
        coordinates_buf[i * 3 + 2] = p_arr[i].pos.z;

        if (flag_vel_) {
          velocities_buf[i * 3    ] = p_arr[i].ori.x;
          velocities_buf[i * 3 + 1] = p_arr[i].ori.y;
          velocities_buf[i * 3 + 2] = p_arr[i].ori.z;
        }

        if (flag_ori_) {
          const double vi_para_module = x_axis_new.dot(p_arr[i].ori);
          const Vec_3<double> vi_perp(p_arr[i].ori - x_axis_new * vi_para_module);
          vel_para_module[i] = vi_para_module;
          theta[i] = std::atan2(z_axis_new.dot(vi_perp), y_axis_new.dot(vi_perp));
        }
      }

      stat = nc_put_vara(ncid_, coordinates_id_, startset, countset, coordinates_buf);
      check_err(stat, __LINE__, __FILE__);
      if (flag_vel_) {
        stat = nc_put_vara(ncid_, velocities_id_, startset, countset, velocities_buf);
        check_err(stat, __LINE__, __FILE__);
      }
      if (flag_ori_) {
        size_t startset1[2] = { time_idx_[0], 0 };
        size_t countset1[2] = { 1, spatial_len_ };
        stat = nc_put_vara(ncid_, vel_mean_id_, startset1, countset1, &vel_mean.x);
        check_err(stat, __LINE__, __FILE__);

        size_t startset2[2] = { time_idx_[0], particle_begin };
        size_t countset2[2] = { 1, particle_count};
        stat = nc_put_vara(ncid_, vel_parallel_id_, startset2, countset2, vel_para_module);
        check_err(stat, __LINE__, __FILE__);
        stat = nc_put_vara(ncid_, theta_id_, startset2, countset2, theta);
        check_err(stat, __LINE__, __FILE__);
      }

      delete[] coordinates_buf;
      delete[] velocities_buf;
      delete[] vel_para_module;
      delete[] theta;
    }
    time_idx_[0] += 1;
    nc_sync(ncid_);
  }
}
