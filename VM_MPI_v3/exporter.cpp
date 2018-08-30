#include "vect.h"
#include "comn.h"
#include "exporter.h" 
#include "mpi.h"
#include "netcdf.h"
#ifndef _MSC_VER
#include "netcdf_par.h"
#endif

Vec_3<double> gl_l;
Vec_3<int> domain_sizes;
double rho0;
int gl_n_par;
double eta;
double eps;
int n_step;
int my_proc;
int tot_proc;
unsigned long long seed;
std::string folder;
std::string base_name;

void ini_output(int gl_np, double eta0, double eps0, int steps, unsigned long long sd,
                const Vec_3<double> &gl_l0, const Vec_3<int>& domain_sizes0) {
  gl_n_par = gl_np;
  eta = eta0;
  eps = eps0;
  n_step = steps;
  seed = sd;
  gl_l = gl_l0;
  domain_sizes = domain_sizes0;
  rho0 = gl_np / (gl_l.x * gl_l.y * gl_l.z);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
  folder = "data" + delimiter;
  if (my_proc == 0) {
    mkdir(folder);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  char str[100];
  snprintf(str, 100, "%g_%.2f_%.3f_%.1f_%llu", gl_l.x, eta, eps, rho0, seed);
  base_name = str;
}

// check whether there is error when outputting netcdf file
void check_err(const int stat, const int line, const char * file) {
  if (stat != NC_NOERR) {
    (void)fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    fflush(stderr);
    exit(1);
  }
}

int get_start_particle_num(int particle_num) {
  int *particle_count_arr = new int[tot_proc];
  MPI_Gather(&particle_num, 1, MPI_INT, particle_count_arr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  int particle_start_num;
  int *particle_start_arr = new int[tot_proc];
  if (my_proc == 0) {
    particle_start_arr[0] = 0;
    for (int i = 1; i < tot_proc; i++) {
      particle_start_arr[i] = particle_start_arr[i - 1] + particle_count_arr[i - 1];
    }
  }
  MPI_Scatter(particle_start_arr, 1, MPI_INT, &particle_start_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
  delete[] particle_count_arr;
  delete[] particle_start_arr;
  return particle_start_num;
}

LogExporter::LogExporter(int interval)
  : BaseLogExporter(folder + base_name + ".log", gl_n_par/tot_proc, n_step, interval){
  fout_ << "\n-------- Parameters --------";
  fout_ << "\nParticle number=" << gl_n_par;
  fout_ << "\nPacking fraction=" << rho0;
  fout_ << "\ndomain lengths=" << gl_l;
  fout_ << "\nblock sizes=" << domain_sizes;
  fout_ << "\nn_step=" << n_step;
  fout_ << "\nseed=" << seed;
  fout_ << "\n\n-------- RUN --------";
  fout_ << "\ntime step\telapsed time" << std::endl;
}

OrderParaExporter::OrderParaExporter(int interval)
  : BaseExporter(n_step, interval) {
  if (my_proc == 0) {
    fout_.open(folder + "phi_" + base_name + ".dat");
  }
  gl_np_ = gl_n_par;
  my_proc_ = my_proc;
}

OrderParaExporter::~OrderParaExporter() {
  if (my_proc_ == 0) {
    fout_.close();
  }
}

FieldExporter::FieldExporter(int frame_interval, int first_frame, int bin_len,
                           const Vec_3<int>& domain_rank,
                           const Vec_3<int>& gl_cells_size,
                           const Vec_3<int>& my_cells_size)
  : frame_len_(NC_UNLIMITED), cg_box_len_(bin_len) {
  time_idx_[0] = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_proc_);
  set_lin_frame(frame_interval, n_step, first_frame);
  char filename[100];
  snprintf(filename, 100, "%sfield_%s.nc", folder.c_str(), base_name.c_str());
#ifdef _MSC_VER
  auto stat = nc_create(filename, NC_NETCDF4, &ncid_);
#else
  auto stat = nc_create_par(filename, NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid_);
#endif
  check_err(stat, __LINE__, __FILE__);

  set_coarse_grain_box(gl_cells_size, my_cells_size, domain_rank);

  int frame_dim;
  int spatial_dim;
  int gl_field_dims[3];

  /* define dimentions */
  stat = nc_def_dim(ncid_, "frame", frame_len_, &frame_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "spatial", spatial_len_, &spatial_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "global_field_z", gl_field_len_[0], &gl_field_dims[0]);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "global_field_y", gl_field_len_[1], &gl_field_dims[1]);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "global_field_x", gl_field_len_[2], &gl_field_dims[2]);
  check_err(stat, __LINE__, __FILE__);

  /* define variables */
  int spatial_dims[1] = { spatial_dim };
  stat = nc_def_var(ncid_, "spatial", NC_CHAR, 1, spatial_dims, &spatial_id_);
  check_err(stat, __LINE__, __FILE__);
  int time_dims[1] = { frame_dim };
  stat = nc_def_var(ncid_, "time", NC_INT, 1, time_dims, &time_id_);
  check_err(stat, __LINE__, __FILE__);
  int den_dims[4] = { frame_dim, gl_field_dims[0], gl_field_dims[1], gl_field_dims[2] };
  stat = nc_def_var(ncid_, "density_field", NC_USHORT, 4, den_dims, &densities_id_);
  check_err(stat, __LINE__, __FILE__);
  int vel_dims[5] = { frame_dim, spatial_dim, gl_field_dims[0], gl_field_dims[1], gl_field_dims[2]};
  stat = nc_def_var(ncid_, "spatial_field", NC_FLOAT, 5, vel_dims, &velocities_id_);
  check_err(stat, __LINE__, __FILE__);

#ifndef _MSC_VER
  stat = nc_var_par_access(ncid_, densities_id_, NC_COLLECTIVE);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_var_par_access(ncid_, velocities_id_, NC_COLLECTIVE);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_var_par_access(ncid_, time_id_, NC_COLLECTIVE);
  check_err(stat, __LINE__, __FILE__);
#endif

  /* assign global attributes */
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "eta", NC_DOUBLE, 1, &eta);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "epsilon", NC_DOUBLE, 1, &eps);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "rho_0", NC_DOUBLE, 1, &rho0);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_ulonglong(ncid_, NC_GLOBAL, "seed", NC_UINT64, 1, &seed);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "Lx", NC_DOUBLE, 1, &gl_l.x);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "Ly", NC_DOUBLE, 1, &gl_l.y);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "Lz", NC_DOUBLE, 1, &gl_l.z);
  check_err(stat, __LINE__, __FILE__);

  /* assign per-variable attributes */
  stat = nc_put_att_int(ncid_, time_id_, "frame_inteval", NC_INT, 1, &frame_interval_);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_int(ncid_, densities_id_, "box_len", NC_INT, 1, &bin_len);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_int(ncid_, velocities_id_, "box_len", NC_INT, 1, &bin_len);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_put_var(ncid_, spatial_id_, "xyz");
  check_err(stat, __LINE__, __FILE__);
}

FieldExporter::~FieldExporter() {
  const auto stat = nc_close(ncid_);
  check_err(stat, __LINE__, __FILE__);
}


void FieldExporter::set_coarse_grain_box(const Vec_3<int>& gl_cells_size,
                                        const Vec_3<int>& my_cells_size,
                                        const Vec_3<int>& domain_rank) {
  if (my_cells_size.x % cg_box_len_ == 0 && 
      my_cells_size.y % cg_box_len_ == 0 && 
      my_cells_size.z % cg_box_len_ == 0) {
    gl_field_len_[2] = gl_cells_size.x / cg_box_len_;
    gl_field_len_[1] = gl_cells_size.y / cg_box_len_;
    gl_field_len_[0] = gl_cells_size.z / cg_box_len_;

    size_t my_field_len[3]{my_cells_size.z / cg_box_len_,
                           my_cells_size.y / cg_box_len_,
                           my_cells_size.x / cg_box_len_};

    den_start_set_[0] = vel_start_set_[0] = 0;
    vel_start_set_[1] = 0;
    den_start_set_[1] = vel_start_set_[2] = my_field_len[0] * domain_rank.z;
    den_start_set_[2] = vel_start_set_[3] = my_field_len[1] * domain_rank.y;
    den_start_set_[3] = vel_start_set_[4] = my_field_len[2] * domain_rank.x;

    den_count_set_[0] = vel_count_set_[0] = 1;
    vel_count_set_[1] = spatial_len_;
    den_count_set_[1] = vel_count_set_[2] = my_field_len[0];
    den_count_set_[2] = vel_count_set_[3] = my_field_len[1];
    den_count_set_[3] = vel_count_set_[4] = my_field_len[2];
    origin_.x = double(den_start_set_[3] * cg_box_len_);
    origin_.y = double(den_start_set_[2] * cg_box_len_);
    origin_.z = double(den_start_set_[1] * cg_box_len_);
  } else {
    std::cerr << "cells size " << gl_cells_size << "are not divisible by box length = "
              << cg_box_len_ << std::endl;
    exit(-1);
  }
}

ParticleExporter::ParticleExporter(int frame_interval, int first_frame, bool flag_vel, bool flag_ori)
                                   : frame_len_(NC_UNLIMITED), atom_len_(gl_n_par),
                                     flag_vel_(flag_vel), flag_ori_(flag_ori) {
  set_lin_frame(frame_interval, n_step, first_frame);
  cell_lengths_data_[0] = gl_l.x;
  cell_lengths_data_[1] = gl_l.y;
  cell_lengths_data_[2] = gl_l.z;
  cell_origin_data_[0] = cell_origin_data_[1] = cell_origin_data_[2] = 0;
  time_idx_[0] = 0;
  gl_np_ = gl_n_par;

  /* open file */
  char filename[100];
  snprintf(filename, 100, "%ssnap_%s.nc", folder.c_str(), base_name.c_str());
#ifdef _MSC_VER
  auto stat = nc_create(filename, NC_NETCDF4, &ncid_);
#else
  auto stat = nc_create_par(filename, NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid_);
#endif
  check_err(stat, __LINE__, __FILE__);

  /* dimension ids */
  int frame_dim;
  int spatial_dim;
  int atom_dim;
  int cell_spatial_dim;

  /* define dimensions */
  stat = nc_def_dim(ncid_, "frame", frame_len_, &frame_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "spatial", spatial_len_, &spatial_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "atom", atom_len_, &atom_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "cell_spatial", cell_spatial_len_, &cell_spatial_dim);
  check_err(stat, __LINE__, __FILE__);

  /* define variables */
  int spatial_dims[1] = { spatial_dim };
  stat = nc_def_var(ncid_, "spatial", NC_CHAR, 1, spatial_dims, &spatial_id_);
  check_err(stat, __LINE__, __FILE__);

  int cell_spatial_dims[1] = { cell_spatial_dim };
  stat = nc_def_var(ncid_, "cell_spatial", NC_CHAR, 1, cell_spatial_dims, &cell_spatial_id_);
  check_err(stat, __LINE__, __FILE__);

  int time_dims[1] = { frame_dim };
  stat = nc_def_var(ncid_, "time", NC_INT, 1, time_dims, &time_id_);
  check_err(stat, __LINE__, __FILE__);

  int cell_lengths_dims[2] = { frame_dim, cell_spatial_dim };
  stat = nc_def_var(ncid_, "cell_lengths", NC_FLOAT, 2, cell_lengths_dims, &cell_lengths_id_);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_var(ncid_, "cell_origin", NC_FLOAT, 2, cell_lengths_dims, &cell_origin_id_);
  check_err(stat, __LINE__, __FILE__);

  int coordinates_dims[3] = { frame_dim, atom_dim, spatial_dim };
  stat = nc_def_var(ncid_, "coordinates", NC_FLOAT, 3, coordinates_dims, &coordinates_id_);
  check_err(stat, __LINE__, __FILE__);

  if (flag_vel_) {
    stat = nc_def_var(ncid_, "velocities", NC_FLOAT, 3, coordinates_dims, &velocities_id_);
    check_err(stat, __LINE__, __FILE__);
  }

  if (flag_ori_) {
    int vel_mean_dims[2] = { frame_dim, spatial_dim };
    stat = nc_def_var(ncid_, "mean_velocities", NC_DOUBLE, 2, vel_mean_dims, &vel_mean_id_);
    check_err(stat, __LINE__, __FILE__);
    int vel_parallel_dims[2] = { frame_dim, atom_dim };
    stat = nc_def_var(ncid_, "v_parallel", NC_FLOAT, 2, vel_parallel_dims, &vel_parallel_id_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var(ncid_, "theta_for_v_perp", NC_FLOAT, 2, vel_parallel_dims, &theta_id_);
    check_err(stat, __LINE__, __FILE__);
  }

#ifndef _MSC_VER
  stat = nc_var_par_access(ncid_, time_id_, NC_COLLECTIVE);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_var_par_access(ncid_, cell_lengths_id_, NC_COLLECTIVE);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_var_par_access(ncid_, cell_origin_id_, NC_COLLECTIVE);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_var_par_access(ncid_, coordinates_id_, NC_COLLECTIVE);
  check_err(stat, __LINE__, __FILE__);
  if (flag_vel_) {
    stat = nc_var_par_access(ncid_, velocities_id_, NC_COLLECTIVE);
    check_err(stat, __LINE__, __FILE__);
  }
  if (flag_ori_) {
    stat = nc_var_par_access(ncid_, vel_mean_id_, NC_COLLECTIVE);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_var_par_access(ncid_, vel_parallel_id_, NC_COLLECTIVE);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_var_par_access(ncid_, theta_id_, NC_COLLECTIVE);
    check_err(stat, __LINE__, __FILE__);
  }
#endif

  /* assign global attributes */
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "title", 18, "Vicsek model in 3d");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "application", 5, "AMBER");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "program", 6, "sander");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "programVersion", 3, "9.0");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "Conventions", 5, "AMBER");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "ConventionVersion", 3, "1.0");
  check_err(stat, __LINE__, __FILE__);

  stat = nc_put_att_double(ncid_, NC_GLOBAL, "eta", NC_DOUBLE, 1, &eta);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "epsilon", NC_DOUBLE, 1, &eps);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "rho_0", NC_DOUBLE, 1, &rho0);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_ulonglong(ncid_, NC_GLOBAL, "seed", NC_UINT64, 1, &seed);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "Lx", NC_DOUBLE, 1, &gl_l.x);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "Ly", NC_DOUBLE, 1, &gl_l.y);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "Lz", NC_DOUBLE, 1, &gl_l.z);
  check_err(stat, __LINE__, __FILE__);

  /* assign per-variable attributes */
  stat = nc_put_att_int(ncid_, time_id_, "frame_inteval", NC_INT, 1, &frame_interval);
  check_err(stat, __LINE__, __FILE__);

  /* leave define mode */
  stat = nc_enddef(ncid_);
  check_err(stat, __LINE__, __FILE__);

  /* assign variable data */
  stat = nc_put_var(ncid_, spatial_id_, "xyz");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_var(ncid_, cell_spatial_id_, "abc");
  check_err(stat, __LINE__, __FILE__);
}

ParticleExporter::~ParticleExporter() {
  const auto stat = nc_close(ncid_);
  check_err(stat, __LINE__, __FILE__);
}

