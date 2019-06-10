#include "domain2D.h"
#include <mpi.h>

Domain_2::Domain_2(const Vec_2<double>& gl_l, const Vec_2<int>& domains_size):
                   gl_l_(gl_l), gl_half_l_(0.5 * gl_l){
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  l_.x = gl_l.x / domains_size.x;
  l_.y = gl_l.y / domains_size.y;
  rank_ = Vec_2<int>(my_rank % domains_size.x, my_rank / domains_size.x);
  origin_.x = l_.x * rank_.x;
  origin_.y = l_.y * rank_.y;
  flag_comm_.x = domains_size.x != 1;
  flag_comm_.y = domains_size.y != 1;
}
