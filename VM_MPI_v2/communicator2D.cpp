#include "communicator2D.h"

void find_shell(const Vec_2<int>& n, const Vec_2<int>& thickness, Vec_2<block_t> shell[2]) {
  for (int ori = 0; ori < 2; ori++) {
    if (thickness[ori]) {
      for (int dim = 0; dim < 2; dim++) {
        if (dim == ori) {
          shell[0][ori].beg[dim] = 0;
          shell[0][ori].end[dim] = 1;
          shell[1][ori].beg[dim] = n[dim] - 1;
          shell[1][ori].end[dim] = n[dim];
        } else if (dim > ori) {
          shell[0][ori].beg[dim] = shell[1][ori].beg[dim] = 0;
          shell[0][ori].end[dim] = shell[1][ori].end[dim] = n[dim];
        } else {
          shell[0][ori].beg[dim] = shell[1][ori].beg[dim] = thickness[dim];
          shell[0][ori].end[dim] = shell[1][ori].end[dim] = n[dim] - thickness[dim];
        }
      }
    }
  }
}

Vec_2<int> decompose_domain(const Vec_2<double> &gl_l) {
  int tot_proc, my_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
  Vec_2<int> domains_size{};
  if (gl_l.x == gl_l.y) {
    int nx = int(std::sqrt(tot_proc));
    while (tot_proc % nx != 0) {
      nx--;
    }
    domains_size = Vec_2<int>(nx, tot_proc / nx);
  } else if (gl_l.x > gl_l.y) {
    const double ratio = gl_l.x / gl_l.y;
    if (my_proc == 0) {
      std::cout << "ratio = " << ratio << std::endl;
    }
    int ny = 1;
    while (ny <= std::sqrt(tot_proc)) {
      if (ratio * ny * ny >= tot_proc && tot_proc % ny == 0) {
        domains_size = Vec_2<int>(tot_proc / ny, ny);
        break;
      }
      ny++;
    }
  } else {
    const double ratio = gl_l.y / gl_l.x;
    int nx = 1;
    while (nx <= std::sqrt(tot_proc)) {
      if (ratio * nx * nx >= tot_proc && tot_proc % nx == 0) {
        domains_size = Vec_2<int>(nx, tot_proc / nx);
        break;
      }
      nx++;
    }
  }
  if (my_proc == 0) {
    std::cout << "decompose the domain (" << gl_l.x << ", " << gl_l.y
      << ") into " << domains_size.x << " by "
      << domains_size.y << " subdomains, each of which has size "
      << gl_l.x / domains_size.x << " by " << gl_l.y / domains_size.y
      << std::endl;
  }
  return domains_size;
}

Communicator::Communicator(const Vec_2<double>& gl_l,
                           double rho0,
                           double amplification,
                           const Vec_2<int>& cells_size,
                           const Vec_2<int> &domains_size) {
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc_);
  const Vec_2<int> domain_rank(my_rank_ % domains_size.x, my_rank_ / domains_size.x);

  find_neighbor(domain_rank, domains_size, flag_comm_);
  set_comm_shell(cells_size);

  Vec_2<double> l(gl_l.x / domains_size.x, gl_l.y / domains_size.y);
  max_buf_size_ = get_max_buf_size(rho0, amplification, l);
  for (int i = 0; i < 4; i++) {
    buf_[i] = new double[max_buf_size_];
    buf_size_[i] = max_buf_size_;
  }

  vacant_pos_.reserve(max_buf_size_);
}

void Communicator::find_neighbor(const Vec_2<int>& domain_rank,
                                 const Vec_2<int>& domain_size,
                                 Vec_2<bool> &flag_comm) {
  const int nx = domain_size.x;
  for (int dim = 0; dim < 2; dim++) {
    if (domain_size[dim] > 1) {
      flag_comm[dim] = true;
      Vec_2<int> prev(domain_rank);
      Vec_2<int> next(domain_rank);
      prev[dim] -= 1;
      next[dim] += 1;
      if (prev[dim] < 0)
        prev[dim] += domain_size[dim];
      if (next[dim] >= domain_size[dim])
        next[dim] -= domain_size[dim];
      neighbor_[dim][0] = prev.x + prev.y * nx;
      neighbor_[dim][1] = next.x + next.y * nx;
    } else {
      flag_comm[dim] = false;
      neighbor_[dim][0] = neighbor_[dim][1] = MPI_PROC_NULL;
    }
  }
}

int Communicator::get_max_buf_size(double rho0, double amplification,
                                    const Vec_2<double> &l) const {
  std::vector<double> area;
  if (flag_comm_.x) {
    area.push_back(l.y);
  }
  if (flag_comm_.y) {
    area.push_back(l.x);
  }

  int max_buf_size;
  if (area.empty()) {
    max_buf_size = 0;

  } else {
    std::sort(area.begin(), area.end(), [](double x, double y) {return x > y; });

    int n0 = int(rho0 * area[0] * amplification);
    max_buf_size = 4 * n0;
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) {
      std::cout << "max area = " << area[0] << std::endl;
      std::cout << "max particle number per communication: " << n0 << " particles" << std::endl;
    }
  }
  return max_buf_size;
}

void Communicator::set_comm_shell(const Vec_2<int>& cells_size) {
  Vec_2<int> thickness{};
  for (int dim = 0; dim < 2; dim++) {
    thickness[dim] = flag_comm_[dim] ? 1 : 0;
  }

  find_shell(cells_size, -thickness, inner_shell_);
  for (int ori = 0; ori < 2; ori++) {
    inner_shell_[0][ori].beg += thickness;
    inner_shell_[0][ori].end += thickness;
    inner_shell_[1][ori].beg += thickness;
    inner_shell_[1][ori].end += thickness;
  }

  const Vec_2<int> extended_cells_size = cells_size + thickness * 2;
  find_shell(extended_cells_size, thickness, outer_shell_);
}
