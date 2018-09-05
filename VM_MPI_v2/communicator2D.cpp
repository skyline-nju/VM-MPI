#include "communicator2D.h"
#define TEST 2

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

void set_comm_block(const Vec_2<int>& cells_size, const Vec_2<bool>& flag_comm,
                    Vec_2<block_t> inner_shell[2], Vec_2<block_t> outer_shell[2]) {
  Vec_2<int> thickness{};
  for (int dim = 0; dim < 2; dim++) {
    thickness[dim] = flag_comm[dim] ? 1 : 0;
  }

  find_shell(cells_size, -thickness, inner_shell);
  for (int ori = 0; ori < 2; ori++) {
    inner_shell[0][ori].beg += thickness;
    inner_shell[0][ori].end += thickness;
    inner_shell[1][ori].beg += thickness;
    inner_shell[1][ori].end += thickness;
  }

  const Vec_2<int> extended_cells_size = cells_size + thickness * 2;
  find_shell(extended_cells_size, thickness, outer_shell);

#if TEST == 2
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0) {
    std::cout << "proc = " << my_rank << std::endl;
    std::cout << "n: " << cells_size << std::endl;
    std::cout << "thickness: " << thickness << std::endl;

    std::cout << "outer shell" << std::endl;
    for (int i = 0; i < 2; i++) {
      for (int dim = 0; dim < 2; dim++) {
        std::cout << outer_shell[i][dim].beg << ";\t" << outer_shell[i][dim].end << "\n";
      }
    }
    std::cout << "inner shell" << std::endl;
    for (int i = 0; i < 2; i++) {
      for (int dim = 0; dim < 2; dim++)
        std::cout << inner_shell[i][dim].beg << ";\t" << inner_shell[i][dim].end << "\n";
    }
    std::cout << std::endl;
    std::cout << "----------" << std::endl;
  }
#endif
}
