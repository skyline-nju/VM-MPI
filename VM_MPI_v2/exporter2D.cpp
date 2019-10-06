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
unsigned long long seed;
std::string folder;
std::string base_name;

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
#ifdef BIRTH_DEATH
  if (gl_l.x == gl_l.y) {
    snprintf(str, 100, "%g_%.3f_%g_%.1f_%llu", gl_l.x, eta, eps * 1e8, rho0, seed);
  } else {
    snprintf(str, 100, "%g_%g_%.3f_%g_%.1f_%llu", gl_l.x, gl_l.y, eta, eps * 1e8, rho0, seed);
  }
#else
  if (gl_l.x == gl_l.y) {
    snprintf(str, 100, "%g_%.2f_%.1f_%llu", gl_l.x, eta, rho0, seed);
  } else {
    snprintf(str, 100, "%g_%g_%.2f_%.1f_%llu", gl_l.x, gl_l.y, eta, rho0, seed);
  }
#endif
#endif
  base_name = str;
}


void exporter_ini(int gl_np, double eta0, double eps0, int steps, unsigned long long sd,
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
  //set_multi_nodes();
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
#ifdef POLAR_ALIGN
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

