#include "run2D.h"

void set_particle_num(int gl_par_num, int& my_par_num, int& my_par_num_max, MPI_Comm group_comm) {
#ifdef USE_MPI
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
  if (my_rank > 0) {
    my_par_num = gl_par_num / tot_proc;
  } else {
    my_par_num = gl_par_num - gl_par_num / tot_proc * (tot_proc - 1);
  }
  my_par_num = my_par_num;
  if (tot_proc > 1) {
    my_par_num_max = 10 * my_par_num;
  } else {
    my_par_num_max = 1.1 * my_par_num;
  }
#else
  my_par_num = gl_par_num;
  my_par_num_max = gl_par_num;
#endif
}


