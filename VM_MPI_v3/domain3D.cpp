#define TEST 10
#include "domain3D.h"

//Domain_3::Domain_3(const Vec3d & gl_l)
//  : l_(gl_l), gl_l_(gl_l), gl_half_l_(gl_l * 0.5) {}
//
//Domain_3::Domain_3(const Vec3d& gl_l, const Vec3i& gl_size,
//                   const Vec3i& gl_cells_size, Vec3i& cells_size,
//                   Vec3d& origin, Vec3b &flag_c)
//  : l_(), gl_l_(gl_l), gl_half_l_(gl_l_ * 0.5), gl_size_(gl_size),
//    gl_cells_size_(gl_cells_size), cells_size_(cells_size) {
//  find_neighbor(rank_, flag_comm_, neighbor);
//  flag_c = flag_comm_;
//  set_l(gl_cells_size, cells_size, l_, origin_);
//  origin = origin_;
//  set_comm_block(cells_size, flag_comm_, inner_shell, outer_shell);
//}
//
//Domain_3::Domain_3(const Vec3d& gl_l, const Vec3i& gl_cells_size,
//                   Vec3i& cells_size, Vec3d& origin, Vec3b& flag_c)
//  : Domain_3(gl_l, partition(gl_l_), gl_cells_size, cells_size, origin, flag_c) {}

void Domain_3::find_neighbor(Vec3i &rank, Vec3b &flag_comm, int neighbor[3][2]) const {
  int my_rank;
  MPI_Comm_rank(group_comm_, &my_rank);
  const int nx = gl_size_.x;
  const int nx_ny = gl_size_.x * gl_size_.y;
  rank.x = (my_rank % nx_ny) % nx;
  rank.y = (my_rank % nx_ny) / nx;
  rank.z = my_rank / nx_ny;
  for (int dim = 0; dim < 3; dim++) {
    if (gl_size_[dim] > 1) {
      flag_comm[dim] = true;
      Vec3i prev(rank_);
      Vec3i next(rank_);
      prev[dim] -= 1;
      next[dim] += 1;
      if (prev[dim] < 0)
        prev[dim] += gl_size_[dim];
      if (next[dim] >= gl_size_[dim])
        next[dim] -= gl_size_[dim];
      neighbor[dim][0] = prev.x + prev.y * nx + prev.z * nx_ny;
      neighbor[dim][1] = next.x + next.y * nx + next.z * nx_ny;
    } else {
      flag_comm[dim] = false;
      neighbor[dim][0] = neighbor[dim][1] = MPI_PROC_NULL;
    }
  }
#if TEST == 0
  if (my_rank == 0) {
    std::cout << "global size: " << gl_size_ << std::endl;
    std::cout << "rank: " << rank_ << std::endl;
    std::cout << "neighbor\n";
    for (int dim = 0; dim < 3; dim ++) {
      std::cout << neighbor[dim][0] << "\t" << neighbor[dim][1] << std::endl;
    }    
  }
#endif
}

void Domain_3::set_l(const Vec3i& gl_cells_size, Vec3i& cells_size,
                     Vec3d& l, Vec3d& origin) const {
  for(int dim = 0; dim < 3; dim++) {
    int pos = 0;
    for (int rank = 0; rank < gl_size_[dim]; rank++) {
      const int nc = (gl_cells_size[dim] - pos) / (gl_size_[dim] - rank);
      if (rank == rank_[dim]) {
        cells_size[dim] = nc;
        break;
      }
      pos += nc;
    }
    const double lc = gl_l_[dim] / gl_cells_size[dim];
    l[dim] = lc * cells_size[dim];
    origin[dim] = lc * pos;
  }
#if TEST == 0
  int my_rank;
  MPI_Comm_rank(group_comm_, &my_rank);
  if (my_rank == 1) {
    std::cout << "origin:\t" << origin << std::endl;
    std::cout << "l:\t" << l << std::endl;
  }
#endif
}

void Domain_3::set_max_buf_size(int gl_par_num, double amplification) {
  std::vector<double> area;
  if (flag_comm_.x) {
    area.push_back(l_.y * l_.z);  
  }
  if (flag_comm_.y) {
    area.push_back(l_.x * l_.z);
  }
  if (flag_comm_.z) {
    area.push_back(l_.x * l_.y);
  }

  if (area.empty()) {
    max_buf_size_ = 0;
    
  } else {
    std::sort(area.begin(), area.end(), [](double x, double y) {return x > y; });
    
    const double rho0 = gl_par_num / (gl_l_.x * gl_l_.y * gl_l_.z);
    int n0 = int(rho0 * area[0] * amplification);
    max_buf_size_ = 6 * n0;
    int my_rank;
    MPI_Comm_rank(group_comm_, &my_rank);

    if (my_rank == 0) {
      std::cout << "max area = " << area[0] << std::endl;
      std::cout << "max particle number per communication: " << n0  << " particles" << std::endl;
    }
  }
}

Vec_3<int> Domain_3::partition(const Vec3d& l, int n_proc) {
  if (n_proc == 1) {
    return Vec_3<int>(1, 1, 1);
  }
  std::vector<int> idx{ 0, 1, 2 };
  std::sort(idx.begin(), idx.end(), [&l](int a, int b) {return l[a] < l[b]; });
  Vec3d l_s{l[idx[0]], l[idx[1]], l[idx[2]]};

  const auto cal_area = [&l_s](int nx, int ny, int nz) {
    double area = 0;
    if (nx > 1) {
      area += l_s.y * l_s.z / (ny * nz);
    }
    if (ny > 1) {
      area += l_s.x * l_s.z / (nx * nz);
    }
    if (nz > 1) {
      area += l_s.x * l_s.y / (nx * ny);
    }
    return area;
  };

  double area_min = -1;
  Vec3i n_opt{};
  for (int i = 1; i < n_proc; i++) {
    if (n_proc % i == 0) {
      int jk = n_proc / i;
      for (int j = i; j < n_proc; j++) {
        if (jk % j == 0) {
          int k = jk / j;
          if (k >= j) {
            const double area = cal_area(i, j, k);
            if (area_min < 0) {
              area_min = area;
              n_opt = Vec3i(i, j, k);
            } else if (area < area_min) {
              area_min = area;
              n_opt = Vec3i(i, j, k);
            }
          }
        }
      }
    }
  }

  return Vec3i(n_opt[idx[0]], n_opt[idx[1]], n_opt[idx[2]]);
}

Vec_3<int> Domain_3::partition(const Vec3d& l, MPI_Comm group_comm) {
#ifdef USE_MPI
  int tot_proc, my_rank;
  MPI_Comm_size(group_comm, &tot_proc);
  MPI_Comm_size(group_comm, &my_rank);
  Vec_3<int> domains_size = partition(l, tot_proc);
  if (my_rank == 0)
    std::cout << "domains size = " << domains_size << std::endl;
  return domains_size;
#else
  return Vec3i(1, 1, 1);
#endif
}

