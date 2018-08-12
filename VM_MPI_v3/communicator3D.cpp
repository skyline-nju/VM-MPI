#include "mpi.h"
#include "communicator3D.h"

Vec_3<int> divide_cubic_domain() {
  int tot_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &tot_rank);
  Vec_3<int> n_domains;
  if (tot_rank % 2 == 1 && tot_rank > 4) {
    std::cout << "tot_rank = " << tot_rank << std::endl;
    exit(1);
  }
  if(tot_rank <= 4) {
    n_domains = Vec_3<int>(1, 1, tot_rank);
  } else {
    switch(tot_rank) {
    case 6:
      n_domains = Vec_3<int>(1, 2, 3);
      break;
    case 8:
      n_domains = Vec_3<int>(2, 2, 2);
      break;
    case 10:
      n_domains = Vec_3<int>(1, 2, 5);
      break;
    case 12:
      n_domains = Vec_3<int>(2, 2, 3);
      break;
    case 24:
      n_domains = Vec_3<int>(2, 3, 4);
      break;
    case 36:
      n_domains = Vec_3<int>(3, 3, 4);
      break;
    case 48:
      n_domains = Vec_3<int>(3, 4, 4);
      break;
    case 60:
      n_domains = Vec_3<int>(3, 4, 5);
      break;
    case 72:
      n_domains = Vec_3<int>(3, 4, 6);
      break;
    case 84:
      n_domains = Vec_3<int>(3, 4, 7);
      break;
    case 96:
      n_domains = Vec_3<int>(4, 4, 6);
      break;
    default:
      std::cout << "tot_rank = " << tot_rank << std::endl;
      exit(1);
    }
  }
  return n_domains;
}

void find_neighbor_proc(const Vec_3<int>& domain_size,
                        Vec_3<int> &domain_rank,
                        int neighbor_proc[3][2],
                        Vec_3<bool>& flag_ext,
                        int test_rank) {
  int tot_rank, my_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &tot_rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

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
    std::cout << domain_rank << std::endl;
    for (int dim = 0; dim < 3; dim++) {
      std::cout << neighbor_proc[dim][0] << "\t" << neighbor_proc[dim][1] << std::endl;
    }
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
  int my_rank = domain_rank.x + domain_rank.y * domain_size.x
    + domain_rank.z * domain_size.x * domain_size.y;
  if (my_rank == 5) {
    std::cout << "domain size: " << domain_size << std::endl;
    std::cout << "domain rank: " << domain_rank << std::endl;
    std::cout << "cells size: " << my_cell_size << std::endl;
    std::cout << "first cell rank: " << my_first_cell << std::endl;
  }
}
