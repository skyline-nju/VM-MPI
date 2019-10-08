#pragma once

#include "config.h"
#include "comn.h"
#include "vect.h"
#include "mpi.h"
#include "cellList2D.h"

void exporter_ini(int gl_np, double eta0, double eps0, int steps, unsigned long long sd,
                const Vec_2<double> &gl_l0, const Vec_2<int> &domain_sizes0);

void create_output_folder(const std::string& fd);

std::string set_base_name();

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




