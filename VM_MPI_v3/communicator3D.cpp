#include "communicator3D.h"

int get_proc_num() {
  int tot_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);
  return tot_proc;
}

int get_proc_rank() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

void find_neighbor_proc(const Vec_3<int>& domain_size,
                        Vec_3<int> &domain_rank,
                        int neighbor_proc[3][2],
                        Vec_3<bool>& flag_ext,
                        int test_rank) {
  const int my_rank = get_proc_rank();
  const int nxny_domains = domain_size.x * domain_size.y;
  domain_rank = Vec_3<int>((my_rank % nxny_domains) % domain_size.x,
                           (my_rank % nxny_domains) / domain_size.x,
                            my_rank / nxny_domains);
  for (int dim = 0; dim < 3; dim++) {
    if (domain_size[dim] > 1) {
      flag_ext[dim] = true;
      int prev = domain_rank[dim] - 1;
      int next = domain_rank[dim] + 1;
      if (prev < 0)
        prev += domain_size[dim];
      if (next >= domain_size[dim])
        next -= domain_size[dim];
      if (dim == 0) {
        neighbor_proc[dim][0] = prev + domain_rank.y * domain_size.x + domain_rank.z * nxny_domains;
        neighbor_proc[dim][1] = next + domain_rank.y * domain_size.x + domain_rank.z * nxny_domains;
      } else if (dim == 1) {
        neighbor_proc[dim][0] = domain_rank.x + prev * domain_size.x + domain_rank.z * nxny_domains;
        neighbor_proc[dim][1] = domain_rank.x + next * domain_size.x + domain_rank.z * nxny_domains;
      } else {
        neighbor_proc[dim][0] = domain_rank.x + domain_rank.y * domain_size.x + prev * nxny_domains;
        neighbor_proc[dim][1] = domain_rank.x + domain_rank.y * domain_size.x + next * nxny_domains;
      }
    } else {
      flag_ext[dim] = false;
      neighbor_proc[dim][0] = MPI_PROC_NULL;
      neighbor_proc[dim][1] = MPI_PROC_NULL;
    }
  }

  if (my_rank == test_rank) {
    std::cout << "domain size:\t" << domain_size << std::endl;
  }
}

void divide_domain(const Vec_3<int>& domain_size,
                   const Vec_3<int>& domain_rank,
                   const Vec_3<int>& gl_cell_size,
                   Vec_3<int>& my_cell_size,
                   Vec_3<int>& my_first_cell) {
  for(int dim = 0; dim < 3; dim++) {
    for (int rank = 0, res = gl_cell_size[dim]; rank < domain_size[dim]; rank++) {
      const int nc = res / (domain_size[dim] - rank);
      if (rank == domain_rank[dim]) {
        my_cell_size[dim] = nc;
        my_first_cell[dim] = gl_cell_size[dim] - res;
      }
      res -= nc;
    }
  }
  const int my_rank = get_proc_rank();
  if (my_rank == 5) {
    std::cout << "domain size: " << domain_size << std::endl;
    std::cout << "domain rank: " << domain_rank << std::endl;
    std::cout << "cells size: " << my_cell_size << std::endl;
    std::cout << "first cell rank: " << my_first_cell << std::endl;
  }
}

int divide_par_evenly(int gl_np) {
  const int tot_proc = get_proc_num();
  const int my_proc = get_proc_rank();
  for (int i = 0, res = gl_np; i < tot_proc; i++) {
    int np = res / (tot_proc - i);
    if (i == my_proc) {
      return np;
    }
    res -= np;
  }
  return 0;
}

Vec_3<block_t> inner_block[2];
Vec_3<block_t> outer_block[2];

void set_comm_block(const Vec_3<int> &cell_size, const Vec_3<bool> &flag_comm) {
  Vec_3<int> n = cell_size;
  int y_beg = flag_comm.y ? 1 : 0;
  int z_beg = flag_comm.z ? 1 : 0;
  int y_end = n.y - y_beg;
  int z_end = n.z - z_beg;
  if (flag_comm.x) {
    inner_block[0].x.beg = Vec_3<int>(1, y_beg, z_beg);
    inner_block[0].x.end = Vec_3<int>(1, y_end, z_end);
    inner_block[1].x.beg = Vec_3<int>(n.x - 2, y_beg, z_beg);
    inner_block[1].x.end = Vec_3<int>(n.x - 2, y_end, z_end);

    outer_block[0].x.beg = Vec_3<int>(0, y_beg, z_beg);
    outer_block[0].x.end = Vec_3<int>(0, y_end, z_end);
    outer_block[1].x.beg = Vec_3<int>(n.x - 1, y_beg, z_beg);
    outer_block[1].x.end = Vec_3<int>(n.x - 1, z_beg, z_end);
  }

  int x_beg = 0;
  int x_end = n.x - x_beg;

  if (flag_comm.y) {
    inner_block[0].y.beg = Vec_3<int>(x_beg, 1, z_beg);
    inner_block[0].y.end = Vec_3<int>(x_end, 1, z_end);
    inner_block[1].y.beg = Vec_3<int>(x_beg, n.y - 2, z_beg);
    inner_block[1].y.end = Vec_3<int>(x_end, n.y - 2, z_end);

    outer_block[0].y.beg = Vec_3<int>(x_beg, 0, z_beg);
    outer_block[0].y.end = Vec_3<int>(x_end, 0, z_end);
    outer_block[1].y.beg = Vec_3<int>(x_beg, n.y - 1, z_beg);
    outer_block[1].y.end = Vec_3<int>(x_end, n.y - 1, z_end);
  }

  y_beg = 0;
  y_end = n.y - y_beg;

  if (flag_comm.z) {
    inner_block[0].z.beg = Vec_3<int>(x_beg, y_beg, 1);
    inner_block[0].z.end = Vec_3<int>(x_end, y_end, 1);
    inner_block[1].z.beg = Vec_3<int>(x_beg, y_beg, n.z - 2);
    inner_block[1].z.end = Vec_3<int>(x_end, y_end, n.z - 2);

    outer_block[0].z.beg = Vec_3<int>(x_beg, y_beg, 0);
    outer_block[0].z.end = Vec_3<int>(x_end, y_end, 0);
    outer_block[1].z.beg = Vec_3<int>(x_beg, y_beg, n.z - 1);
    outer_block[1].z.end = Vec_3<int>(x_end, y_end, n.z - 1);
  }
}

const Vec_3<block_t>* get_inner_block() {
  return inner_block;
}

const Vec_3<block_t>* get_outer_block() {
  return outer_block;
}

