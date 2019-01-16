#define TEST 10
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
//Domain_2::Domain_2(const Vec2d & gl_l)
//  : l_(gl_l), gl_l_(gl_l), gl_half_l_(gl_l * 0.5) {}
//
//Domain_2::Domain_2(const Vec2d& gl_l, const Vec2i& gl_size,
//                   const Vec2i& gl_cells_size, Vec2i& cells_size,
//                   Vec2d& origin, Vec2b &flag_c)
//                   : l_(), gl_l_(gl_l), gl_half_l_(gl_l_ * 0.5),
//                     gl_size_(gl_size),
//  gl_cells_size_(gl_cells_size), cells_size_(cells_size) {
//  find_neighbor(rank_, flag_comm_, neighbor);
//  flag_c = flag_comm_;
//  set_l(gl_cells_size, cells_size, l_, origin_);
//  origin = origin_;
//  set_comm_block(cells_size, flag_comm_, inner_shell, outer_shell);
//}
//
//Domain_2::Domain_2(const Vec2d& gl_l, const Vec2i& gl_cells_size,
//  Vec2i& cells_size, Vec2d& origin, Vec2b& flag_c)
//  : Domain_2(gl_l, partition(gl_l_), gl_cells_size, cells_size, origin, flag_c) {
//}

//void Domain_2::find_neighbor(Vec2i &rank, Vec2b &flag_comm, int neighbor[2][2]) const {
//  int my_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//  const int nx = gl_size_.x;
//  rank.x = my_rank % nx;
//  rank.y = my_rank / nx;
//  for (int dim = 0; dim < 2; dim++) {
//    if (gl_size_[dim] > 1) {
//      flag_comm[dim] = true;
//      Vec2i prev(rank_);
//      Vec2i next(rank_);
//      prev[dim] -= 1;
//      next[dim] += 1;
//      if (prev[dim] < 0)
//        prev[dim] += gl_size_[dim];
//      if (next[dim] >= gl_size_[dim])
//        next[dim] -= gl_size_[dim];
//      neighbor[dim][0] = prev.x + prev.y * nx;
//      neighbor[dim][1] = next.x + next.y * nx;
//    } else {
//      flag_comm[dim] = false;
//      neighbor[dim][0] = neighbor[dim][1] = MPI_PROC_NULL;
//    }
//  }
//}

//void Domain_2::set_l(const Vec2i& gl_cells_size, Vec2i& cells_size,
//                     Vec2d& l, Vec2d& origin) const {
//  for (int dim = 0; dim < 2; dim++) {
//    int pos = 0;
//    for (int rank = 0; rank < gl_size_[dim]; rank++) {
//      const int nc = (gl_cells_size[dim] - pos) / (gl_size_[dim] - rank);
//      if (rank == rank_[dim]) {
//        cells_size[dim] = nc;
//        break;
//      }
//      pos += nc;
//    }
//    const double lc = gl_l_[dim] / gl_cells_size[dim];
//    l[dim] = lc * cells_size[dim];
//    origin[dim] = lc * pos;
//  }
//}

//void Domain_2::set_max_buf_size(int gl_par_num, double amplification) {
//  std::vector<double> area;
//  if (flag_comm_.x) {
//    area.push_back(l_.y);
//  }
//  if (flag_comm_.y) {
//    area.push_back(l_.x);
//  }
//
//  if (area.empty()) {
//    max_buf_size_ = 0;
//
//  } else {
//    std::sort(area.begin(), area.end(), [](double x, double y) {return x > y; });
//
//    const double rho0 = gl_par_num / (gl_l_.x * gl_l_.y);
//    int n0 = int(rho0 * area[0] * amplification);
//    max_buf_size_ = 4 * n0;
//    int my_rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//
//    if (my_rank == 0) {
//      std::cout << "max area = " << area[0] << std::endl;
//      std::cout << "max particle number per communication: " << n0 << " particles" << std::endl;
//    }
//  }
//}



//Vec_2<int> Domain_2::partition(const Vec2d& l, int n_proc) {
//  Vec_2<int> domain_size{};
//  if (n_proc == 1) {
//    domain_size = Vec_2<int>(1, 1);
//  } else if (l.x > l.y) {
//    double ratio = l.x / l.y;
//    if (ratio >= n_proc) {
//      domain_size = Vec_2<int>(n_proc, 1);
//    } else if (ratio * 2 >= n_proc) {
//      domain_size = Vec_2<int>(n_proc / 2, 2);
//    }
//  } else if (l.y / l.x >= n_proc) {
//    domain_size = Vec_2<int>(1, n_proc);
//  } else {
//    int nx = int(std::sqrt(n_proc));
//    while (n_proc % nx != 0) {
//      nx -= 1;
//    }
//    int ny = n_proc / nx;
//    domain_size = Vec_2<int>(nx, ny);
//  }
//  return domain_size;
//}

//Vec_2<int> Domain_2::partition(const Vec2d& l) {
//#ifdef USE_MPI
//  int tot_proc, my_rank;
//  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
//  MPI_Comm_size(MPI_COMM_WORLD, &my_rank);
//  Vec_2<int> domains_size = partition(l, tot_proc);
//  if (my_rank == 0)
//    std::cout << "domains size = " << domains_size << std::endl;
//  return domains_size;
//#else
//  return Vec2i(1, 1);
//#endif
//}

