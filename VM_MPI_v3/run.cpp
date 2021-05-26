#include "run.h"

int ini_my_par_num(int gl_par_num, MPI_Comm group_comm) {
#ifdef USE_MPI
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
  int my_par_num;
  if (my_rank > 0) {
    my_par_num = gl_par_num / tot_proc;
  } else {
    my_par_num = gl_par_num - gl_par_num / tot_proc * (tot_proc - 1);
  }
  return my_par_num;
#else
  return gl_par_num;
#endif
}
