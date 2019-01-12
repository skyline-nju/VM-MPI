#include "exporter2D.h"

Vec_2<double> gl_l;
Vec_2<int> domain_sizes;
double rho0;
int gl_n_par;
double eta;
double eps;
int n_step;
int my_proc;
int tot_proc;
int t_half = 0;
int t_eq;
double H0;
unsigned long long seed;
std::string folder;
std::string base_name;

int my_host = 0;
#ifdef NP_PER_NODE
int tot_host;
MPI_Group gl_group;
MPI_Group *host_group = nullptr;
MPI_Comm *host_comm = nullptr;
#endif

void create_output_folder() {
  folder = "data" + delimiter;
  if (my_proc == 0) {
    mkdir(folder);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void set_base_name() {
  char str[100];
#ifdef DISORDER_ON
  if (gl_l.x == gl_l.y) {
    snprintf(str, 100, "%g_%.2f_%.3f_%.1f_%llu", gl_l.x, eta, eps, rho0, seed);
  } else {
    snprintf(str, 100, "%g_%g_%.2f_%.3f_%.1f_%llu", gl_l.x, gl_l.y, eta, eps, rho0, seed);
  }
#else
  if (gl_l.x == gl_l.y) {
    snprintf(str, 100, "%g_%.2f_%.1f_%llu", gl_l.x, eta, rho0, seed);
  } else {
    snprintf(str, 100, "%g_%g_%.2f_%.1f_%llu", gl_l.x, gl_l.y, eta, rho0, seed);
  }
#endif
  base_name = str;
}

void set_mult_nodes() {
  if (tot_proc > NP_PER_NODE) {
    tot_host = tot_proc / NP_PER_NODE;
    host_group = new MPI_Group[tot_host];
    host_comm = new MPI_Comm[tot_host];
    MPI_Comm_group(MPI_COMM_WORLD, &gl_group);
    for (int i = 0; i < tot_host; i++) {
      int *rank_arr = new int[NP_PER_NODE];
      for (int j = 0; j < NP_PER_NODE; j++) {
        rank_arr[j] = j + i * NP_PER_NODE;
      }
      MPI_Group_incl(gl_group, NP_PER_NODE, rank_arr, &host_group[i]);
      delete[] rank_arr;
    }

    for (int i = 0; i < tot_host; i++) {
      MPI_Comm_create(MPI_COMM_WORLD, host_group[i], &host_comm[i]);
    }
    my_host = my_proc / NP_PER_NODE;
  }
}

void ini_output(int gl_np, double eta0, double eps0, int steps, unsigned long long sd,
                const Vec_2<double> &gl_l0, const Vec_2<int> &domain_sizes0) {
  gl_n_par = gl_np;
  eta = eta0;
  eps = eps0;
  n_step = steps;
  seed = sd;
  gl_l = gl_l0;
  domain_sizes = domain_sizes0;
  rho0 = gl_np / (gl_l.x * gl_l.y);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);

  create_output_folder();
  set_base_name();
  set_mult_nodes();
}

void output_finalize() {
  if (tot_proc > NP_PER_NODE) {
    for (int i = 0; i < tot_host; i++) {
      if(MPI_GROUP_NULL != host_group[i]) {
        MPI_Group_free(&host_group[i]);
      }
      if (MPI_COMM_NULL != host_comm[i]) {
        MPI_Comm_free(&host_comm[i]);
      }
    }
    delete[] host_group;
    delete[] host_comm;
  }
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
  : BaseLogExporter(folder + base_name + ".log", gl_n_par / tot_proc, n_step, interval) {
  fout_ << "\n-------- Parameters --------";
  fout_ << "\nParticle number=" << gl_n_par;
  fout_ << "\nPacking fraction=" << rho0;
  fout_ << "\ndomain lengths=" << gl_l;
  fout_ << "\nblock sizes=" << domain_sizes;
  fout_ << "\nn_step=" << n_step;
  fout_ << "\nseed=" << seed;

#ifdef SCALAR_NOISE
  fout_ << "\nscalar noise";
#else
  fout_ << "\nvectorial noise";
#endif
  fout_ << " = " << eta;

#ifdef DISORDER_ON
#ifdef RANDOM_TORQUE
  fout_ << "random torque";
#elif RANDOM_FIELD
  fout_ << "random filed";
#else
  fout_ << "random stress";
#endif
  fout_ << "  = " << eps;
#endif
#ifdef POLAR_ALIGNE
  fout_ << "\npolar alignment";
#else
  fout_ << "\nnematic alignment";
#endif
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
                             const Vec_2<int>& domain_rank,
                             const Vec_2<int>& gl_cells_size,
                             const Vec_2<int>& my_cells_size)
                             : frame_len_(NC_UNLIMITED), cg_box_len_(bin_len) {
  time_idx_[0] = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_proc_);
  set_lin_frame(frame_interval, n_step, first_frame);
  snprintf(filename_, 100, "%sfield_%s_host%d.nc", folder.c_str(), base_name.c_str(), my_host);

  int stat;
#ifdef _MSC_VER
  stat = nc_create(filename_, NC_NETCDF4, &ncid_);
#else
  if (tot_proc <= NP_PER_NODE) {
    stat = nc_create_par(filename_, NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid_);
  } else {
    stat = nc_create_par(filename_, NC_NETCDF4 | NC_MPIIO, host_comm[my_host], MPI_INFO_NULL, &ncid_);
  }
#endif

  check_err(stat, __LINE__, __FILE__);

  set_coarse_grain_box(gl_cells_size, my_cells_size, domain_rank);

  int frame_dim;
  int spatial_dim;
  int gl_field_dims[2];

  /* define dimentions */
  stat = nc_def_dim(ncid_, "frame", frame_len_, &frame_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "spatial", spatial_len_, &spatial_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "NY", gl_field_len_[0], &gl_field_dims[0]);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "NX", gl_field_len_[1], &gl_field_dims[1]);
  check_err(stat, __LINE__, __FILE__);

  /* define variables */
  int spatial_dims[1] = { spatial_dim };
  stat = nc_def_var(ncid_, "spatial", NC_CHAR, 1, spatial_dims, &spatial_id_);
  check_err(stat, __LINE__, __FILE__);
  int time_dims[1] = { frame_dim };
  stat = nc_def_var(ncid_, "time", NC_INT, 1, time_dims, &time_id_);
  check_err(stat, __LINE__, __FILE__);
  int den_dims[3] = { frame_dim, gl_field_dims[0], gl_field_dims[1]};
  stat = nc_def_var(ncid_, "density_field", NC_USHORT, 3, den_dims, &densities_id_);
  check_err(stat, __LINE__, __FILE__);
  int vel_dims[4] = { frame_dim, spatial_dim, gl_field_dims[0], gl_field_dims[1]};
  stat = nc_def_var(ncid_, "velocity_field", NC_FLOAT, 4, vel_dims, &velocities_id_);
  check_err(stat, __LINE__, __FILE__);

  int n_host_id;
  int host_rank_id;
  if (tot_proc > NP_PER_NODE) {
    stat = nc_def_var(ncid_, "host_size", NC_INT, 1, spatial_dims, &n_host_id);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var(ncid_, "host_rank", NC_INT, 1, spatial_dims, &host_rank_id);
    check_err(stat, __LINE__, __FILE__);
  }

  /* assign global attributes */
  if (t_half > 0) {
    stat = nc_put_att_text(ncid_, NC_GLOBAL, "title", 38, "Vicsek model in 2d");
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_int(ncid_, NC_GLOBAL, "t_half", NC_INT, 1, &t_half);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_int(ncid_, NC_GLOBAL, "t_equil", NC_INT, 1, &t_eq);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "h0", NC_DOUBLE, 1, &H0);
    check_err(stat, __LINE__, __FILE__);
  } else {
    stat = nc_put_att_text(ncid_, NC_GLOBAL, "title", 38, "Vicsek model in 2d");
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "epsilon", NC_DOUBLE, 1, &eps);
    check_err(stat, __LINE__, __FILE__);
  }
#ifdef SCALAR_NOISE
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "noise", 5, "scalar");
#else
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "noise", 9, "vectorial");
#endif
  check_err(stat, __LINE__, __FILE__);
#ifdef POLAR_ALIGN
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "alignment", 5, "polar");
#else
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "alignment", 7, "nematic");
#endif
  check_err(stat, __LINE__, __FILE__);


#ifdef RANDOM_TORQUE
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "disorder", 2, "RT");
  check_err(stat, __LINE__, __FILE__);
#elif RANDOM_FIELD
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "disorder", 2, "RF");
  check_err(stat, __LINE__, __FILE__);
#elif RANDOM_STRESS
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "disorder", 2, "RS");
  check_err(stat, __LINE__, __FILE__);
#endif
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "eta", NC_DOUBLE, 1, &eta);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "rho_0", NC_DOUBLE, 1, &rho0);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_ulonglong(ncid_, NC_GLOBAL, "seed", NC_UINT64, 1, &seed);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "Lx", NC_DOUBLE, 1, &gl_l.x);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "Ly", NC_DOUBLE, 1, &gl_l.y);
  check_err(stat, __LINE__, __FILE__);

  /* assign per-variable attributes */
  stat = nc_put_att_int(ncid_, time_id_, "frame_inteval", NC_INT, 1, &frame_interval_);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_int(ncid_, densities_id_, "box_len", NC_INT, 1, &bin_len);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_int(ncid_, velocities_id_, "box_len", NC_INT, 1, &bin_len);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_put_var(ncid_, spatial_id_, "xy");
  check_err(stat, __LINE__, __FILE__);

  if (tot_proc > NP_PER_NODE) {
    int n_host[2] = { n_host_.y, n_host_.x };
    stat = nc_put_var(ncid_, n_host_id, n_host);
    check_err(stat, __LINE__, __FILE__);
    int host_rank[2] = { my_host / n_host_.x, my_host % n_host_.x };
    stat = nc_put_var(ncid_, host_rank_id, host_rank);
    check_err(stat, __LINE__, __FILE__);
  }
  stat = nc_close(ncid_);
  check_err(stat, __LINE__, __FILE__);
}

void FieldExporter::open() {
  int stat;
#ifdef _MSC_VER
  stat = nc_open(filename_, NC_WRITE, &ncid_);
#else
  if (tot_proc <= NP_PER_NODE) {
  stat = nc_open_par(filename_, NC_WRITE | NC_MPIIO, MPI_COMM_WORLD,
                     MPI_INFO_NULL, &ncid_);
  } else {
  stat = nc_open_par(filename_, NC_WRITE | NC_MPIIO, host_comm[my_host],
                          MPI_INFO_NULL, &ncid_);
  }
#endif
  check_err(stat, __LINE__, __FILE__);
}

void FieldExporter::set_coarse_grain_box(const Vec_2<int>& gl_cells_size,
                                         const Vec_2<int>& my_cells_size,
                                         const Vec_2<int>& domain_rank) {
  if (my_cells_size.x % cg_box_len_ == 0 &&
      my_cells_size.y % cg_box_len_ == 0) {
    n_host_.x = n_host_.y = 1;
    if (tot_proc > NP_PER_NODE) {
      if (domain_sizes.x <= NP_PER_NODE && NP_PER_NODE % domain_sizes.x == 0) {
          n_host_.x = 1;
      } else if (domain_sizes.x > NP_PER_NODE && domain_sizes.x % NP_PER_NODE == 0) {
          n_host_.x = domain_sizes.x / NP_PER_NODE;   
      } else {
        if (my_proc == 0) {
          std::cerr << "NP_PER_NODE = " << NP_PER_NODE
                    << "domain_sizes.x = " << domain_sizes.x << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
      n_host_.y = tot_host / n_host_.x;
    }
    gl_field_len_[1] = gl_cells_size.x / cg_box_len_ / n_host_.x;
    gl_field_len_[0] = gl_cells_size.y / cg_box_len_ / n_host_.y;

    size_t my_field_len[2]{my_cells_size.y / cg_box_len_,
                           my_cells_size.x / cg_box_len_ };

    den_start_set_[0] = vel_start_set_[0] = 0;
    vel_start_set_[1] = 0;
    den_start_set_[1] = vel_start_set_[2] = my_field_len[0] * domain_rank.y % gl_field_len_[0];
    den_start_set_[2] = vel_start_set_[3] = my_field_len[1] * domain_rank.x % gl_field_len_[1];

    den_count_set_[0] = vel_count_set_[0] = 1;
    vel_count_set_[1] = spatial_len_;
    den_count_set_[1] = vel_count_set_[2] = my_field_len[0];
    den_count_set_[2] = vel_count_set_[3] = my_field_len[1];
    origin_.x = double(my_field_len[1] * domain_rank.x * cg_box_len_);
    origin_.y = double(my_field_len[0] * domain_rank.y * cg_box_len_);

  } else {
    if (my_proc == 0) {
      std::cerr << "cells size " << my_cells_size << " are not divisible by box length = "
                << cg_box_len_ << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
}