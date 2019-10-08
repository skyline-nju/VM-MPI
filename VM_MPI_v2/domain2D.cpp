#include "domain2D.h"
#include <mpi.h>

Grid_2::Grid_2(const Vec_2<double>& gl_l, const Vec_2<int>& proc_size, 
               const Vec_2<int>& proc_rank, double r_cut) {
  gl_n_.x = int(gl_l.x / r_cut);
  gl_n_.y = int(gl_l.y / r_cut);

  lc_.x = gl_l.x / gl_n_.x;
  lc_.y = gl_l.y / gl_n_.y;

  Vec_2<int> n_per_proc(gl_n_.x / proc_size.x, gl_n_.y / proc_size.y);
 
  origin_ = proc_rank * n_per_proc;
  if (proc_rank.x == proc_size.x - 1) {
    n_.x = gl_n_.x - proc_rank.x * n_per_proc.x;
  } else {
    n_.x = n_per_proc.x;
  }
  if (proc_rank.y == proc_size.y - 1) {
    n_.y = gl_n_.y - proc_rank.y * n_per_proc.y;
  } else {
    n_.y = n_per_proc.y;
  }
}

SubDomain_2::SubDomain_2(const Vec_2<double>& gl_l, const Vec_2<int>& proc_size, double r_cut) :
  size_(proc_size), rank_(get_proc_rank(size_)), grid_(gl_l, size_, rank_, r_cut),
  flag_comm_(size_.x != 1, size_.y != 1), gl_l_(gl_l), half_gl_l_(0.5 * gl_l_) {

  l_.x = grid_.n().x * grid_.lc().x;
  l_.y = grid_.n().y * grid_.lc().y;

  origin_.x = grid_.origin().x * grid_.lc().x;
  origin_.y = grid_.origin().y * grid_.lc().y;
}

void SubDomain_2::find_neighbor(int(*neighbor)[2]) const {
  const int nx = size_.x;
  for (int dim = 0; dim < 2; dim++) {
    if (size_[dim] > 1) {
      Vec_2<int> prev(rank_);
      Vec_2<int> next(rank_);
      prev[dim] -= 1;
      next[dim] += 1;
      if (prev[dim] < 0)
        prev[dim] += size_[dim];
      if (next[dim] >= size_[dim])
        next[dim] -= size_[dim];
      neighbor[dim][0] = prev.x + prev.y * nx;
      neighbor[dim][1] = next.x + next.y * nx;
    } else {
      neighbor[dim][0] = neighbor[dim][1] = MPI_PROC_NULL;
    }
  }
}

const Vec_2<int> SubDomain_2::decompose(const Vec_2<double>& gl_l) {
  int tot_proc, my_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
  Vec_2<int> proc_size{};
  if (gl_l.x == gl_l.y) {
    int nx = int(std::sqrt(tot_proc));
    while (tot_proc % nx != 0) {
      nx--;
    }
    proc_size = Vec_2<int>(nx, tot_proc / nx);
  } else if (gl_l.x > gl_l.y) {
    const double ratio = gl_l.x / gl_l.y;
    if (my_proc == 0) {
      std::cout << "ratio = " << ratio << std::endl;
    }
    int ny = 1;
    while (ny <= std::sqrt(tot_proc)) {
      if (ratio * ny * ny >= tot_proc && tot_proc % ny == 0) {
        proc_size = Vec_2<int>(tot_proc / ny, ny);
        break;
      }
      ny++;
    }
  } else {
    const double ratio = gl_l.y / gl_l.x;
    int nx = 1;
    while (nx <= std::sqrt(tot_proc)) {
      if (ratio * nx * nx >= tot_proc && tot_proc % nx == 0) {
        proc_size = Vec_2<int>(nx, tot_proc / nx);
        break;
      }
      nx++;
    }
  }
  if (my_proc == 0) {
    std::cout << "decompose the domain (" << gl_l.x << ", " << gl_l.y
      << ") into " << proc_size.x << " by "
      << proc_size.y << " subdomains, each of which has size "
      << gl_l.x / proc_size.x << " by " << gl_l.y / proc_size.y
      << std::endl;
  }
  return proc_size;
}

Vec_2<int> get_proc_rank(const Vec_2<int>& proc_size) {
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  return Vec_2<int>(my_rank % proc_size.x, my_rank / proc_size.x);
}
