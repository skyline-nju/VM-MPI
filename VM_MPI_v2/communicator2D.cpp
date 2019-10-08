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

int Communicator::get_max_buf_size(double rho0, double amplification,
                                   const Vec_2<double>& l) const {
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